import { EXAMPLES, normalizeDNA, parseFASTA, validateOptions } from './core.js';
const $ = (id) => document.getElementById(id),
  el = (tag, text) => {
    const n = document.createElement(tag);
    if (text !== undefined) n.textContent = text;
    return n;
  },
  canvas = $('matrix'),
  ctx = canvas.getContext('2d');
let worker = null,
  revision = 0,
  importRevision = 0,
  result = null,
  names = ['Reference', 'Query'],
  pixels = null,
  cellRow = 0,
  cellColumn = 0;
function counts() {
  for (const id of ['reference', 'query']) {
    try {
      $(id + '-count').textContent = `${normalizeDNA($(id).value).length} / 500 bases`;
    } catch (error) {
      $(id + '-count').textContent = error.message;
    }
  }
}
function invalidate(message = 'Sequences or scoring changed. Align again to update the result.') {
  revision++;
  worker?.terminate();
  worker = null;
  result = null;
  pixels = null;
  $('output').hidden = true;
  $('cancel').disabled = true;
  $('status').textContent = message;
  counts();
}
function read() {
  const options = { mode: $('mode').value, reverse: $('reverse').checked };
  for (const key of ['match', 'mismatch', 'gapOpen', 'gapExtend']) {
    if ($(key).value.trim() === '') throw Error('Fill every scoring parameter.');
    options[key] = Number($(key).value);
  }
  return {
    reference: normalizeDNA($('reference').value),
    query: normalizeDNA($('query').value),
    options: validateOptions(options),
  };
}
$('form').addEventListener('input', (event) => {
  if (event.target.type === 'file') return;
  importRevision++;
  invalidate();
});
$('form').onsubmit = (event) => {
  event.preventDefault();
  run();
};
$('cancel').onclick = () => invalidate('Alignment cancelled.');
function run() {
  let data;
  try {
    data = read();
  } catch (error) {
    invalidate(error.message);
    return;
  }
  invalidate('Aligning in a background worker…');
  const ticket = revision;
  $('cancel').disabled = false;
  try {
    worker = new Worker('./worker.js', { type: 'module' });
    worker.onmessage = ({ data }) => {
      if (ticket !== revision) return;
      worker.terminate();
      worker = null;
      $('cancel').disabled = true;
      if (data.error) {
        invalidate(data.error);
        return;
      }
      result = data.result;
      display(data.milliseconds);
    };
    worker.onerror = (event) => invalidate(`Alignment failed: ${event.message}`);
    worker.postMessage(data);
  } catch (error) {
    invalidate(error.message);
  }
}
function sequenceLine(label, start, text, end) {
  const line = el('div');
  line.append(`${label} ${String(start).padStart(4, ' ')}  `);
  for (const base of text) {
    const span = el('span', base);
    span.className = base === '-' ? 'base-gap' : `base-${base.toLowerCase()}`;
    line.append(span);
  }
  line.append(`  ${end}`);
  return line;
}
function display(milliseconds) {
  const r = result;
  $('output').hidden = false;
  const editorHadFocus = $('sequence-editor').contains(document.activeElement);
  $('sequence-editor').open = false;
  if (editorHadFocus) $('sequence-editor').querySelector('summary').focus({ preventScroll: true });
  $('score').textContent = r.score;
  $('identity').textContent = r.identity === null ? '—' : `${(100 * r.identity).toFixed(1)}%`;
  $('gaps').textContent = `${r.gaps} / ${r.gapRuns}`;
  $('ambiguous').textContent = r.ambiguous;
  $('cigar').textContent = r.cigar || '(empty alignment)';
  $('status').textContent =
    `${r.options.mode === 'global' ? 'Global' : 'Local'} alignment complete · ${milliseconds.toFixed(2)} ms · traceback score verified.`;
  $('ranges').textContent =
    `Reference [${r.referenceRange.join(', ')}) · Query in ${r.options.reverse ? 'reverse-complement' : 'forward'} orientation [${r.queryRange.join(', ')}) · Original query coordinates [${r.originalQueryRange.join(', ')}). Identity excludes gaps and ambiguous pairs.`;
  $('alignment').replaceChildren();
  let ref = r.referenceRange[0],
    query = r.queryRange[0];
  for (let offset = 0; offset < r.alignedReference.length; offset += 40) {
    const a = r.alignedReference.slice(offset, offset + 40),
      b = r.alignedQuery.slice(offset, offset + 40),
      refEnd = ref + a.replaceAll('-', '').length,
      queryEnd = query + b.replaceAll('-', '').length,
      block = el('pre');
    block.append(sequenceLine('REF', ref, a, refEnd));
    const symbols = [...a]
      .map((base, i) =>
        base === '-' || b[i] === '-'
          ? ' '
          : base === 'N' || b[i] === 'N'
            ? ':'
            : base === b[i]
              ? '|'
              : '.',
      )
      .join('');
    block.append(el('div', '          ' + symbols));
    block.append(sequenceLine('QRY', query, b, queryEnd));
    $('alignment').append(block);
    ref = refEnd;
    query = queryEnd;
  }
  if (!r.alignedReference.length)
    $('alignment').append(
      el(
        'p',
        r.options.mode === 'local'
          ? 'No positively scoring local alignment.'
          : 'Both sequences are empty.',
      ),
    );
  canvas.width = r.columns;
  canvas.height = r.rows;
  pixels = ctx.createImageData(r.columns, r.rows);
  let low = Infinity,
    high = -Infinity;
  for (const value of r.matrix) {
    low = Math.min(low, value);
    high = Math.max(high, value);
  }
  for (let i = 0; i < r.matrix.length; i++) {
    const t = (r.matrix[i] - low) / (high - low || 1),
      base = i * 4;
    pixels.data[base] = Math.round(239 - 174 * t);
    pixels.data[base + 1] = Math.round(236 - 129 * t);
    pixels.data[base + 2] = Math.round(243 - 82 * t);
    pixels.data[base + 3] = 255;
  }
  $('row').max = r.rows - 1;
  $('column').max = r.columns - 1;
  cellRow = r.referenceRange[1];
  cellColumn = r.queryRange[1];
  paintMatrix();
}
function paintMatrix() {
  if (!result || !pixels) return;
  ctx.putImageData(pixels, 0, 0);
  ctx.beginPath();
  result.trace.forEach((p, index) => {
    index ? ctx.lineTo(p.column + 0.5, p.row + 0.5) : ctx.moveTo(p.column + 0.5, p.row + 0.5);
  });
  ctx.strokeStyle = '#ffffff';
  ctx.lineWidth = Math.max(0.35, result.columns / 350);
  ctx.lineJoin = 'round';
  ctx.stroke();
  ctx.strokeStyle = '#e59079';
  ctx.lineWidth = Math.max(0.2, result.columns / 700);
  ctx.strokeRect(cellColumn + 0.08, cellRow + 0.08, 0.84, 0.84);
  $('row').value = cellRow;
  $('column').value = cellColumn;
  $('cell').textContent =
    `Prefix pair (${cellRow}, ${cellColumn}) · best-state score ${result.matrix[cellRow * result.columns + cellColumn]}`;
}
function pointer(event) {
  if (!result) return;
  const rect = canvas.getBoundingClientRect();
  cellColumn = Math.max(
    0,
    Math.min(
      result.columns - 1,
      Math.floor(((event.clientX - rect.left) / rect.width) * result.columns),
    ),
  );
  cellRow = Math.max(
    0,
    Math.min(result.rows - 1, Math.floor(((event.clientY - rect.top) / rect.height) * result.rows)),
  );
  paintMatrix();
}
canvas.onpointermove = pointer;
canvas.onpointerdown = (event) => {
  canvas.focus();
  pointer(event);
};
canvas.onkeydown = (event) => {
  if (!result) return;
  const moves = { ArrowLeft: [0, -1], ArrowRight: [0, 1], ArrowUp: [-1, 0], ArrowDown: [1, 0] };
  if (moves[event.key]) {
    event.preventDefault();
    const [di, dj] = moves[event.key];
    cellRow = Math.max(0, Math.min(result.rows - 1, cellRow + di));
    cellColumn = Math.max(0, Math.min(result.columns - 1, cellColumn + dj));
    paintMatrix();
  }
};
$('inspect').onclick = () => {
  if (!result) return;
  const row = Number($('row').value),
    column = Number($('column').value);
  if (
    !$('row').value ||
    !$('column').value ||
    !Number.isInteger(row) ||
    !Number.isInteger(column) ||
    row < 0 ||
    column < 0 ||
    row >= result.rows ||
    column >= result.columns
  ) {
    $('cell').textContent = 'Choose valid integer prefix lengths.';
    return;
  }
  cellRow = row;
  cellColumn = column;
  paintMatrix();
};
function loadExample() {
  importRevision++;
  const key = $('example').value,
    data = EXAMPLES[key];
  $('reference').value = data.reference;
  $('query').value = data.query;
  $('mode').value = key === 'local' ? 'local' : 'global';
  $('reverse').checked = key === 'reverse';
  names = ['Reference', 'Query'];
  $('reference-label').textContent = 'REFERENCE / A';
  $('query-label').textContent = 'QUERY / B';
  counts();
  run();
}
$('example').onchange = loadExample;
$('fasta').onchange = async (event) => {
  const file = event.target.files[0];
  event.target.value = '';
  if (!file) return;
  const ticket = ++importRevision;
  try {
    if (file.size > 20000) throw Error('FASTA file must be smaller than 20 KB.');
    const records = parseFASTA(await file.text());
    if (ticket !== importRevision) return;
    $('reference').value = records[0].sequence;
    $('query').value = records[1].sequence;
    names = records.map((r) => r.name);
    $('reference-label').textContent = `REFERENCE / ${names[0]}`;
    $('query-label').textContent = `QUERY / ${names[1]}`;
    invalidate('FASTA loaded. Align to calculate the result.');
  } catch (error) {
    if (ticket === importRevision)
      $('status').textContent = `FASTA import failed: ${error.message}`;
  }
};
function download(data, name, type = 'application/json') {
  const url = URL.createObjectURL(
      new Blob([type === 'application/json' ? JSON.stringify(data, null, 2) : data], { type }),
    ),
    a = el('a');
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}
$('report').onclick = () => {
  if (!result) return;
  const { matrix, ...report } = result;
  download({ ...report, names }, 'sequence-alignment.json');
};
$('aligned-fasta').onclick = () => {
  if (!result) return;
  const wrap = (text) => text.match(/.{1,70}/g)?.join('\n') || '';
  download(
    `>reference start=${result.referenceRange[0]} end=${result.referenceRange[1]}\n${wrap(result.alignedReference)}\n>query orientation=${result.options.reverse ? 'reverse-complement' : 'forward'} start=${result.queryRange[0]} end=${result.queryRange[1]}\n${wrap(result.alignedQuery)}\n`,
    'sequence-alignment.fasta',
    'text/plain',
  );
};
loadExample();
