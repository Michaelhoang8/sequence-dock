export function normalizeDNA(value) {
  if (typeof value !== 'string' || value.length > 10000)
    throw Error('Provide a DNA sequence as text.');
  const sequence = value.replace(/\s/g, '').toUpperCase();
  if (sequence.length > 500 || /[^ACGTN]/.test(sequence))
    throw Error('Sequences must contain at most 500 A, C, G, T, or N bases.');
  return sequence;
}
export function reverseComplement(sequence) {
  return normalizeDNA(sequence)
    .split('')
    .reverse()
    .map((base) => ({ A: 'T', T: 'A', C: 'G', G: 'C', N: 'N' })[base])
    .join('');
}
export function parseFASTA(text) {
  if (typeof text !== 'string' || text.length > 20000)
    throw Error('FASTA input exceeds 20,000 characters.');
  const records = [];
  let current = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith(';')) continue;
    if (line.startsWith('>')) {
      if (records.length >= 2) throw Error('Load exactly two FASTA records.');
      const name = line.slice(1).trim();
      if (!name || name.length > 120) throw Error('FASTA headers need 1–120 characters.');
      current = { name, sequence: '' };
      records.push(current);
    } else {
      if (!current) throw Error('Each FASTA sequence needs a preceding >header.');
      current.sequence += line;
    }
  }
  if (records.length !== 2) throw Error('Load exactly two FASTA records.');
  for (const record of records) {
    record.sequence = normalizeDNA(record.sequence);
    if (!record.sequence) throw Error('FASTA records cannot be empty.');
  }
  return records;
}
export function validateOptions(input = {}) {
  const defaults = {
    mode: 'global',
    match: 2,
    mismatch: 2,
    gapOpen: 3,
    gapExtend: 1,
    reverse: false,
  };
  if (
    !input ||
    typeof input !== 'object' ||
    Array.isArray(input) ||
    Object.keys(input).some((k) => !Object.hasOwn(defaults, k))
  )
    throw Error('Unknown alignment option.');
  const options = { ...defaults, ...input };
  if (!['global', 'local'].includes(options.mode) || typeof options.reverse !== 'boolean')
    throw Error('Invalid alignment mode or orientation.');
  for (const [key, min, max] of [
    ['match', 1, 5],
    ['mismatch', 0, 10],
    ['gapOpen', 1, 10],
    ['gapExtend', 0, 5],
  ])
    if (!Number.isInteger(options[key]) || options[key] < min || options[key] > max)
      throw Error(`${key} must be an integer from ${min} to ${max}.`);
  return options;
}
function substitution(a, b, options) {
  return a === 'N' || b === 'N' ? 0 : a === b ? options.match : -options.mismatch;
}
export function scoreAlignment(a, b, input = {}) {
  const options = validateOptions(input);
  if (
    typeof a !== 'string' ||
    typeof b !== 'string' ||
    a.length !== b.length ||
    /[^ACGTN-]/.test(a + b)
  )
    throw Error('Aligned sequences need equal lengths and DNA or gap symbols.');
  let score = 0,
    previous = '';
  for (let i = 0; i < a.length; i++) {
    if (a[i] === '-' && b[i] === '-') throw Error('A column cannot contain two gaps.');
    const gap = a[i] === '-' ? 'I' : b[i] === '-' ? 'D' : '';
    if (gap) score -= gap === previous ? options.gapExtend : options.gapOpen;
    else score += substitution(a[i], b[i], options);
    previous = gap;
  }
  return score;
}
export function align(reference, query, input = {}) {
  const options = validateOptions(input),
    a = normalizeDNA(reference),
    originalQuery = normalizeDNA(query),
    b = options.reverse ? reverseComplement(originalQuery) : originalQuery;
  const rows = a.length + 1,
    columns = b.length + 1,
    size = rows * columns,
    local = options.mode === 'local',
    scores = Array.from({ length: 3 }, () => new Float64Array(size).fill(local ? 0 : -Infinity)),
    previous = Array.from({ length: 3 }, () => new Int8Array(size).fill(-1));
  scores[0][0] = 0;
  if (!local) {
    for (let i = 1; i < rows; i++) {
      scores[1][i * columns] = -options.gapOpen - (i - 1) * options.gapExtend;
      previous[1][i * columns] = i === 1 ? 0 : 1;
    }
    for (let j = 1; j < columns; j++) {
      scores[2][j] = -options.gapOpen - (j - 1) * options.gapExtend;
      previous[2][j] = j === 1 ? 0 : 2;
    }
  }
  let best = 0,
    endI = 0,
    endJ = 0,
    endState = 0;
  function choose(values, state, at) {
    let value = values[0],
      from = 0;
    for (let k = 1; k < 3; k++)
      if (values[k] > value) {
        value = values[k];
        from = k;
      }
    if (local && value <= 0) {
      value = 0;
      from = -1;
    }
    scores[state][at] = value;
    previous[state][at] = from;
    return value;
  }
  for (let i = 1; i < rows; i++)
    for (let j = 1; j < columns; j++) {
      const at = i * columns + j,
        diag = at - columns - 1,
        up = at - columns,
        left = at - 1,
        s = substitution(a[i - 1], b[j - 1], options);
      choose(
        scores.map((values) => values[diag] + s),
        0,
        at,
      );
      choose(
        [
          scores[0][up] - options.gapOpen,
          scores[1][up] - options.gapExtend,
          scores[2][up] - options.gapOpen,
        ],
        1,
        at,
      );
      choose(
        [
          scores[0][left] - options.gapOpen,
          scores[1][left] - options.gapOpen,
          scores[2][left] - options.gapExtend,
        ],
        2,
        at,
      );
      if (local)
        for (let state = 0; state < 3; state++)
          if (scores[state][at] > best) {
            best = scores[state][at];
            endI = i;
            endJ = j;
            endState = state;
          }
    }
  if (!local) {
    endI = a.length;
    endJ = b.length;
    const at = size - 1;
    endState = 0;
    for (let state = 1; state < 3; state++)
      if (scores[state][at] > scores[endState][at]) endState = state;
    best = scores[endState][at];
  }
  let i = endI,
    j = endJ,
    state = endState;
  const aa = [],
    bb = [],
    trace = [{ row: i, column: j }];
  while (i > 0 || j > 0) {
    const at = i * columns + j;
    if (local && scores[state][at] <= 0) break;
    const prior = previous[state][at];
    if (prior < 0) throw Error('Unexpected traceback boundary.');
    if (state === 0) {
      aa.push(a[--i]);
      bb.push(b[--j]);
    } else if (state === 1) {
      aa.push(a[--i]);
      bb.push('-');
    } else {
      aa.push('-');
      bb.push(b[--j]);
    }
    state = prior;
    trace.push({ row: i, column: j });
  }
  const alignedReference = aa.reverse().join(''),
    alignedQuery = bb.reverse().join('');
  trace.reverse();
  if (scoreAlignment(alignedReference, alignedQuery, options) !== best)
    throw Error('Traceback score failed independent verification.');
  const operations = [];
  let exact = 0,
    mismatches = 0,
    ambiguous = 0,
    gaps = 0,
    gapRuns = 0,
    lastGap = '';
  for (let k = 0; k < aa.length; k++) {
    let op;
    if (aa[k] === '-') {
      op = 'I';
      gaps++;
    } else if (bb[k] === '-') {
      op = 'D';
      gaps++;
    } else if (aa[k] === 'N' || bb[k] === 'N') {
      op = 'M';
      ambiguous++;
    } else if (aa[k] === bb[k]) {
      op = '=';
      exact++;
    } else {
      op = 'X';
      mismatches++;
    }
    if (op === 'I' || op === 'D') {
      if (lastGap !== op) gapRuns++;
      lastGap = op;
    } else lastGap = '';
    if (operations.at(-1)?.op === op) operations.at(-1).count++;
    else operations.push({ op, count: 1 });
  }
  const matrix = new Float64Array(size);
  for (let at = 0; at < size; at++)
    matrix[at] = Math.max(scores[0][at], scores[1][at], scores[2][at]);
  return {
    reference: a,
    query: b,
    originalQuery,
    options,
    score: best,
    alignedReference,
    alignedQuery,
    referenceRange: [i, endI],
    queryRange: [j, endJ],
    originalQueryRange: options.reverse ? [b.length - endJ, b.length - j] : [j, endJ],
    cigar: operations.map((o) => `${o.count}${o.op}`).join(''),
    exact,
    mismatches,
    ambiguous,
    gaps,
    gapRuns,
    identity: exact + mismatches ? exact / (exact + mismatches) : null,
    trace,
    matrix,
    rows,
    columns,
  };
}
export const EXAMPLES = {
  insertion: {
    reference: 'ACGTTGACCTGATCGTAGCTAGGCTAACCGT',
    query: 'ACGTTGACCTGAAATCGTAGCTAGGCTAACCGT',
  },
  local: { reference: 'TTTTACGTCGATGCTAGCAAA', query: 'GGGACGTCGATGCTAGCCCC' },
  reverse: { reference: 'AACCGTATGCTAGG', query: 'CCTAGCATACGGTT' },
};
