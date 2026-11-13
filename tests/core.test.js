import test from 'node:test';
import assert from 'node:assert/strict';
import {
  align,
  scoreAlignment,
  normalizeDNA,
  reverseComplement,
  parseFASTA,
  validateOptions,
} from '../core.js';
function brute(a, b, o, i = 0, j = 0, last = '') {
  if (i === a.length && j === b.length) return 0;
  let best = -Infinity;
  if (i < a.length && j < b.length)
    best =
      (a[i] === 'N' || b[j] === 'N' ? 0 : a[i] === b[j] ? o.match : -o.mismatch) +
      brute(a, b, o, i + 1, j + 1, '');
  if (i < a.length)
    best = Math.max(
      best,
      -(last === 'D' ? o.gapExtend : o.gapOpen) + brute(a, b, o, i + 1, j, 'D'),
    );
  if (j < b.length)
    best = Math.max(
      best,
      -(last === 'I' ? o.gapExtend : o.gapOpen) + brute(a, b, o, i, j + 1, 'I'),
    );
  return best;
}
test('exact global match and empty boundaries have correct scores', () => {
  assert.equal(align('ACGT', 'ACGT').score, 8);
  assert.equal(align('', 'ACGT').score, -6);
  assert.equal(align('ACGT', '').score, -6);
  assert.equal(align('', '').score, 0);
  assert.equal(align('ACGT', 'ACGT').cigar, '4=');
});
test('one gap run costs open plus length minus one extensions', () => {
  const result = align('AAAACCCC', 'AAAAGGGCCCC', {
    match: 3,
    mismatch: 5,
    gapOpen: 4,
    gapExtend: 1,
  });
  assert.equal(result.score, 18);
  assert.equal(result.gapRuns, 1);
  assert.equal(result.gaps, 3);
  assert.equal(result.cigar, '4=3I4=');
});
test('global scores agree with enumeration of every alignment on small inputs', () => {
  let seed = 993;
  const random = (n) => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return seed % n;
  };
  for (let trial = 0; trial < 220; trial++) {
    const a = Array.from({ length: random(5) }, () => 'ACGN'[random(4)]).join(''),
      b = Array.from({ length: random(5) }, () => 'ACGN'[random(4)]).join(''),
      o = {
        match: 1 + random(4),
        mismatch: random(7),
        gapOpen: 1 + random(4),
        gapExtend: random(3),
      };
    const r = align(a, b, o);
    assert.equal(r.score, brute(a, b, o));
    assert.equal(scoreAlignment(r.alignedReference, r.alignedQuery, o), r.score);
    assert.equal(r.alignedReference.replaceAll('-', ''), a);
    assert.equal(r.alignedQuery.replaceAll('-', ''), b);
  }
});
test('local score equals the best global alignment of all substring pairs', () => {
  const options = { mode: 'local', match: 2, mismatch: 3, gapOpen: 2, gapExtend: 1 };
  for (const a of ['ACG', 'GAA', 'NNG', 'TTC'])
    for (const b of ['CGA', 'AAA', 'NGT', 'GCC']) {
      let expected = 0;
      for (let i = 0; i < a.length; i++)
        for (let z = i + 1; z <= a.length; z++)
          for (let j = 0; j < b.length; j++)
            for (let q = j + 1; q <= b.length; q++)
              expected = Math.max(expected, brute(a.slice(i, z), b.slice(j, q), options));
      const r = align(a, b, options);
      assert.equal(r.score, expected);
      assert.equal(r.alignedReference.replaceAll('-', ''), a.slice(...r.referenceRange));
      assert.equal(r.alignedQuery.replaceAll('-', ''), b.slice(...r.queryRange));
    }
});
test('local alignment excludes unrelated ends and no positive match returns empty', () => {
  const r = align('TTACGTAA', 'GGACGTCC', { mode: 'local' });
  assert.equal(r.score, 8);
  assert.deepEqual(r.referenceRange, [2, 6]);
  assert.deepEqual(r.queryRange, [2, 6]);
  assert.equal(r.cigar, '4=');
  const none = align('AAAA', 'CCCC', { mode: 'local' });
  assert.equal(none.score, 0);
  assert.equal(none.alignedReference, '');
});
test('adjacent opposing gaps are considered when cheaper than a mismatch', () => {
  const r = align('A', 'C', { mismatch: 10, gapOpen: 1, gapExtend: 1 });
  assert.equal(r.score, -2);
  assert.equal(r.gapRuns, 2);
});
test('ambiguous N is neutral and excluded from known-base identity', () => {
  const r = align('AN', 'AT');
  assert.equal(r.score, 2);
  assert.equal(r.ambiguous, 1);
  assert.equal(r.identity, 1);
  assert.equal(r.cigar, '1=1M');
  assert.equal(align('NN', 'NN').identity, null);
});
test('reverse complement and original query coordinates are consistent', () => {
  assert.equal(reverseComplement('a cgtn'), 'NACGT');
  const r = align('ACGT', 'TTACGTGG', { mode: 'local', reverse: true });
  assert.equal(r.query, 'CCACGTAA');
  assert.deepEqual(r.queryRange, [2, 6]);
  assert.deepEqual(r.originalQueryRange, [2, 6]);
  assert.equal(r.score, 8);
  const asymmetric = align('ACG', 'TTTCGTGG', { mode: 'local', reverse: true });
  assert.deepEqual(asymmetric.queryRange, [2, 5]);
  assert.deepEqual(asymmetric.originalQueryRange, [3, 6]);
});
test('FASTA parser normalizes two multiline records and rejects malformed inputs', () => {
  assert.deepEqual(parseFASTA('>one\nac gt\n>two\nAC\nGT'), [
    { name: 'one', sequence: 'ACGT' },
    { name: 'two', sequence: 'ACGT' },
  ]);
  for (const text of [
    'ACGT',
    '>a\nACGT',
    '>a\n> b\nACGT',
    '>a\nAC\n>b\nGT\n>c\nAA',
    '>a\nAU\n>b\nGT',
  ])
    assert.throws(() => parseFASTA(text));
});
test('traceback matrix is bounded, deterministic, and follows legal steps', () => {
  const r = align('ACGT', 'AGT'),
    again = align('ACGT', 'AGT');
  assert.deepEqual(r, again);
  assert.equal(r.matrix.length, r.rows * r.columns);
  for (let k = 1; k < r.trace.length; k++) {
    const di = r.trace[k].row - r.trace[k - 1].row,
      dj = r.trace[k].column - r.trace[k - 1].column;
    assert.ok((di === 1 && dj === 1) || (di === 1 && dj === 0) || (di === 0 && dj === 1));
  }
});
test('strict input and scoring validation rejects unsupported alphabets and parameters', () => {
  for (const value of ['A'.repeat(501), 'ACGU', 'AC-G', null])
    assert.throws(() => normalizeDNA(value));
  for (const o of [
    { mode: 'overlap' },
    { gapOpen: 0 },
    { gapExtend: -1 },
    { match: 1.5 },
    { reverse: 1 },
    { extra: 1 },
  ])
    assert.throws(() => validateOptions(o));
  assert.throws(() => scoreAlignment('-', '-'));
});
