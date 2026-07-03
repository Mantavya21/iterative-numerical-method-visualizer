'use strict';

/*FORMATTING HELPERS*/
function fmt(n, digits) {
  digits = digits === undefined ? 6 : digits;
  if (n === null || n === undefined || Number.isNaN(n)) return '\u2014';
  if (!isFinite(n)) return n > 0 ? '\u221e' : '\u2212\u221e';
  if (n !== 0 && (Math.abs(n) < 1e-4 || Math.abs(n) >= 1e6)) return n.toExponential(4);
  var fixed = Number(n.toFixed(digits));
  return fixed.toString();
}

function isComplexObj(v) {
  return !!v && typeof v === 'object' && 're' in v && 'im' in v;
}

function fmtComplex(c, digits) {
  if (c === null || c === undefined) return '\u2014';
  if (typeof c === 'number') return fmt(c, digits);
  if (!isComplexObj(c)) return '\u2014';
  if (Math.abs(c.im) < 1e-9) return fmt(c.re, digits);
  var sign = c.im >= 0 ? '+' : '\u2212';
  return fmt(c.re, digits) + ' ' + sign + ' ' + fmt(Math.abs(c.im), digits) + 'i';
}

function reOf(v) {
  if (v === null || v === undefined) return null;
  if (typeof v === 'number') return v;
  if (isComplexObj(v)) return v.re;
  return null;
}

/*MATH ENGINE*/
function compileFunction(exprStr) {
  try {
    var node = math.parse(exprStr);
    var code = node.compile();
    var f = function (x) { return code.evaluate({ x: x }); };
    var fPrime = null;
    try {
      var dnode = math.derivative(node, 'x');
      var dcode = dnode.compile();
      fPrime = function (x) { return dcode.evaluate({ x: x }); };
    } catch (e) {
      fPrime = null;
    }
    return { f: f, fPrime: fPrime, error: null };
  } catch (e) {
    return { f: null, fPrime: null, error: e.message };
  }
}

function safeEval(fn, x) {
  if (!fn) return null;
  try {
    var v = fn(x);
    if (typeof v !== 'number' || !isFinite(v)) return null;
    return v;
  } catch (e) {
    return null;
  }
}

/*ALGORITHM: NEWTON-RAPHSON*/
function runNewton(exprStr, x0, tol, maxIter) {
  var parsed = compileFunction(exprStr);
  if (parsed.error) return { steps: [], error: 'Could not parse f(x): ' + parsed.error };
  if (!parsed.fPrime) return { steps: [], error: "Could not compute a symbolic derivative for this expression." };
  var f = parsed.f, fPrime = parsed.fPrime;
  var steps = [];
  var x = x0;
  for (var n = 0; n <= maxIter; n++) {
    var fx = safeEval(f, x);
    if (fx === null) {
      steps.push({ n: n, x: x, fx: null, failed: true, errorMsg: 'f(x) is undefined at x = ' + fmt(x) + '.' });
      break;
    }
    var step = { n: n, x: x, fx: fx };
    if (n > 0) {
      var prevX = steps[n - 1].x;
      step.relError = Math.abs(x) > 1e-14 ? Math.abs((x - prevX) / x) * 100 : Math.abs(x - prevX) * 100;
    }
    if (Math.abs(fx) < tol) {
      step.converged = true;
      steps.push(step);
      break;
    }
    var fpx = safeEval(fPrime, x);
    step.fpx = fpx;
    if (fpx === null || Math.abs(fpx) < 1e-10) {
      step.failed = true;
      step.errorMsg = "The derivative f'(x) \u2248 0 at x = " + fmt(x) + " \u2014 the tangent line is (near) horizontal and never meets the x-axis, so Newton\u2013Raphson breaks down here.";
      steps.push(step);
      break;
    }
    var xNext = x - fx / fpx;
    step.xNext = xNext;
    steps.push(step);
    if (!isFinite(xNext)) {
      steps.push({ n: n + 1, x: xNext, fx: null, failed: true, errorMsg: 'The iteration diverged to infinity.' });
      break;
    }
    if (n === maxIter) step.maxedOut = true;
    x = xNext;
  }
  return { steps: steps, error: null };
}

/*ALGORITHM: SECANT*/
function runSecant(exprStr, x0, x1, tol, maxIter) {
  var parsed = compileFunction(exprStr);
  if (parsed.error) return { steps: [], error: 'Could not parse f(x): ' + parsed.error };
  var f = parsed.f;
  var steps = [];
  var f0 = safeEval(f, x0);
  if (f0 === null) {
    steps.push({ n: 0, x: x0, fx: null, failed: true, errorMsg: 'f(x) is undefined at x\u2080 = ' + fmt(x0) + '.' });
    return { steps: steps, error: null };
  }
  steps.push({ n: 0, x: x0, fx: f0 });
  var f1 = safeEval(f, x1);
  if (f1 === null) {
    steps.push({ n: 1, x: x1, fx: null, failed: true, errorMsg: 'f(x) is undefined at x\u2081 = ' + fmt(x1) + '.' });
    return { steps: steps, error: null };
  }
  var rel1 = Math.abs(x1) > 1e-14 ? Math.abs((x1 - x0) / x1) * 100 : Math.abs(x1 - x0) * 100;
  steps.push({ n: 1, x: x1, fx: f1, relError: rel1 });

  var xPrev = x0, fPrev = f0, xCurr = x1, fCurr = f1;
  for (var n = 2; n <= maxIter + 1; n++) {
    if (Math.abs(fCurr) < tol) {
      steps[steps.length - 1].converged = true;
      break;
    }
    var denom = fCurr - fPrev;
    if (Math.abs(denom) < 1e-12) {
      steps.push({ n: n, x: NaN, fx: null, failed: true, errorMsg: 'f(x\u2099) \u2212 f(x\u2099\u208b\u2081) \u2248 0 \u2014 the secant line is (near) horizontal and never meets the x-axis, so the method breaks down here.' });
      break;
    }
    var xNext = xCurr - fCurr * (xCurr - xPrev) / denom;
    var fNext = safeEval(f, xNext);
    if (fNext === null) {
      steps.push({ n: n, x: xNext, fx: null, failed: true, errorMsg: 'f(x) is undefined at x = ' + fmt(xNext) + '.' });
      break;
    }
    var rel = Math.abs(xNext) > 1e-14 ? Math.abs((xNext - xCurr) / xNext) * 100 : Math.abs(xNext - xCurr) * 100;
    var st = { n: n, x: xNext, fx: fNext, relError: rel };
    if (n === maxIter + 1) st.maxedOut = true;
    steps.push(st);
    xPrev = xCurr; fPrev = fCurr; xCurr = xNext; fCurr = fNext;
  }
  return { steps: steps, error: null };
}

/*ALGORITHM: MULLER (complex-capable)*/
function runMuller(exprStr, x0, x1, x2, tol, maxIter) {
  var parsed = compileFunction(exprStr);
  if (parsed.error) return { steps: [], error: 'Could not parse f(x): ' + parsed.error };
  var f = parsed.f;
  function cf(cx) {
    try {
      var v = f(cx);
      return math.complex(v);
    } catch (e) {
      return null;
    }
  }
  var xs = [math.complex(x0, 0), math.complex(x1, 0), math.complex(x2, 0)];
  var steps = [];
  for (var n = 0; n <= maxIter; n++) {
    var xa = xs[0], xb = xs[1], xc = xs[2];
    var fa = cf(xa), fb = cf(xb), fc = cf(xc);
    if (!fa || !fb || !fc) {
      steps.push({ n: n, x: xc, fx: null, failed: true, errorMsg: 'f(x) is undefined at one of the current interpolation points.' });
      break;
    }
    var h0 = math.subtract(xb, xa);
    var h1 = math.subtract(xc, xb);
    if (math.abs(h0) < 1e-14 || math.abs(h1) < 1e-14) {
      steps.push({ n: n, x: xc, fx: fc, failed: true, errorMsg: 'Two of the interpolation points have converged to the same value \u2014 the parabola is undefined.' });
      break;
    }
    var delta0 = math.divide(math.subtract(fb, fa), h0);
    var delta1 = math.divide(math.subtract(fc, fb), h1);
    var a = math.divide(math.subtract(delta1, delta0), math.add(h1, h0));
    var b = math.add(math.multiply(a, h1), delta1);
    var c = fc;
    var disc = math.sqrt(math.subtract(math.multiply(b, b), math.multiply(4, math.multiply(a, c))));
    var denomPlus = math.add(b, disc);
    var denomMinus = math.subtract(b, disc);
    var denom = math.abs(denomPlus) > math.abs(denomMinus) ? denomPlus : denomMinus;
    if (math.abs(denom) < 1e-14) {
      steps.push({ n: n, x: xc, fx: fc, failed: true, errorMsg: 'Both candidate denominators of the Muller update vanish \u2014 the method breaks down here.' });
      break;
    }
    var xNext = math.subtract(xc, math.divide(math.multiply(2, c), denom));
    var fNext = cf(xNext);
    var relError = math.abs(xNext) > 1e-14
      ? math.abs(math.divide(math.subtract(xNext, xc), xNext)) * 100
      : math.abs(math.subtract(xNext, xc)) * 100;
    var step = { n: n, x: xNext, fx: fNext, relError: relError, quadPoints: [xa, xb, xc], quadFs: [fa, fb, fc] };
    if (fNext && math.abs(fNext) < tol) step.converged = true;
    if (n === maxIter) step.maxedOut = true;
    steps.push(step);
    if (step.converged) break;
    xs = [xb, xc, xNext];
  }
  return { steps: steps, error: null };
}

/*ALGORITHM: NATURAL CUBIC SPLINE*/
function runSpline(points) {
  var n = points.length - 1;
  var x = points.map(function (p) { return p.x; });
  var a = points.map(function (p) { return p.y; });
  var h = [];
  for (var i = 0; i < n; i++) h.push(x[i + 1] - x[i]);
  for (i = 0; i < n; i++) {
    if (h[i] <= 0) return { steps: [], error: 'Data points must be sorted with strictly increasing x values.' };
  }
  var alpha = new Array(n + 1).fill(0);
  for (i = 1; i < n; i++) {
    alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1]);
  }
  var l = new Array(n + 1).fill(0);
  var mu = new Array(n + 1).fill(0);
  var z = new Array(n + 1).fill(0);
  l[0] = 1; mu[0] = 0; z[0] = 0;
  for (i = 1; i < n; i++) {
    l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
  }
  l[n] = 1; z[n] = 0;
  var c = new Array(n + 1).fill(0);
  var b = new Array(n).fill(0);
  var d = new Array(n).fill(0);
  for (var j = n - 1; j >= 0; j--) {
    c[j] = z[j] - mu[j] * c[j + 1];
    b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
    d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
  }
  var steps = [];
  for (j = 0; j < n; j++) {
    steps.push({ n: j, x0: x[j], x1: x[j + 1], a: a[j], b: b[j], c: c[j], d: d[j] });
  }
  return { steps: steps, error: null };
}

function evalSplineAt(steps, xq) {
  for (var i = 0; i < steps.length; i++) {
    var s = steps[i];
    if (xq >= s.x0 - 1e-9 && xq <= s.x1 + 1e-9) {
      var dx = xq - s.x0;
      return s.a + s.b * dx + s.c * dx * dx + s.d * dx * dx * dx;
    }
  }
  return null;
}

/*STATE*/
var state = {
  method: 'newton',
  steps: [],
  computeError: null,
  currentIndex: 0,
  playing: false,
  playTimer: null,
  animT: 1,
  raf: null,
  points: [
    { x: 0, y: 1 }, { x: 1, y: 3 }, { x: 2, y: 2 }, { x: 3, y: 5 }, { x: 4, y: 4 }
  ],
  params: {
    newton: { x0: 2, tol: 0.000001, maxIter: 20 },
    secant: { x0: 1, x1: 2, tol: 0.000001, maxIter: 20 },
    muller: { x0: 0, x1: 1, x2: 2, tol: 0.000001, maxIter: 20 },
    spline: { evalX: '' }
  }
};

var paramDefs = {
  newton: [
    { id: 'x0', label: 'Initial guess x\u2080', step: 'any' },
    { id: 'tol', label: 'Tolerance |f(x)| <', step: 'any' },
    { id: 'maxIter', label: 'Max iterations', step: '1' }
  ],
  secant: [
    { id: 'x0', label: 'x\u2080', step: 'any' },
    { id: 'x1', label: 'x\u2081', step: 'any' },
    { id: 'tol', label: 'Tolerance |f(x)| <', step: 'any' },
    { id: 'maxIter', label: 'Max iterations', step: '1' }
  ],
  muller: [
    { id: 'x0', label: 'x\u2080', step: 'any' },
    { id: 'x1', label: 'x\u2081', step: 'any' },
    { id: 'x2', label: 'x\u2082', step: 'any' },
    { id: 'tol', label: 'Tolerance |\u0394x| <', step: 'any' },
    { id: 'maxIter', label: 'Max iterations', step: '1' }
  ],
  spline: [
    { id: 'evalX', label: 'Evaluate at x (optional)', step: 'any', optional: true }
  ]
};

var methodNames = {
  newton: 'NEWTON\u2013RAPHSON',
  secant: 'SECANT',
  muller: "MULLER\u2019S METHOD",
  spline: 'CUBIC SPLINE'
};

/*DOM REFS*/
var dom = {};
function cacheDom() {
  dom.methodTabs = document.getElementById('methodTabs');
  dom.equationGroup = document.getElementById('equationGroup');
  dom.equationInput = document.getElementById('equationInput');
  dom.equationStatus = document.getElementById('equationStatus');
  dom.splinePointsGroup = document.getElementById('splinePointsGroup');
  dom.pointsTable = document.getElementById('pointsTable');
  dom.addPointBtn = document.getElementById('addPointBtn');
  dom.paramRow = document.getElementById('paramRow');
  dom.computeBtn = document.getElementById('computeBtn');
  dom.configNote = document.getElementById('configNote');
  dom.stepCurrent = document.getElementById('stepCurrent');
  dom.stepTotal = document.getElementById('stepTotal');
  dom.canvas = document.getElementById('plotCanvas');
  dom.canvasEmpty = document.getElementById('canvasEmpty');
  dom.errorCard = document.getElementById('errorCard');
  dom.errorCardTitle = document.getElementById('errorCardTitle');
  dom.errorCardBody = document.getElementById('errorCardBody');
  dom.btnReset = document.getElementById('btnReset');
  dom.btnPrev = document.getElementById('btnPrev');
  dom.btnPlay = document.getElementById('btnPlay');
  dom.btnNext = document.getElementById('btnNext');
  dom.playIcon = document.getElementById('playIcon');
  dom.stepSlider = document.getElementById('stepSlider');
  dom.transportLegend = document.getElementById('transportLegend');
  dom.dataTable = document.getElementById('dataTable');
  dom.dataTableHead = document.getElementById('dataTableHead');
  dom.dataTableBody = document.getElementById('dataTableBody');
  dom.tableEmpty = document.getElementById('tableEmpty');
  dom.readoutMethod = document.getElementById('readoutMethod');
  dom.readoutState = document.getElementById('readoutState');
  dom.readoutRoot = document.getElementById('readoutRoot');
}

/*CONFIG RENDERING*/
function renderMethodTabs() {
  var tabs = dom.methodTabs.querySelectorAll('.method-tab');
  tabs.forEach(function (tab) {
    var active = tab.getAttribute('data-method') === state.method;
    tab.classList.toggle('is-active', active);
    tab.setAttribute('aria-selected', active ? 'true' : 'false');
  });
}

function renderParamRow() {
  var defs = paramDefs[state.method];
  var pvals = state.params[state.method];
  dom.paramRow.innerHTML = '';
  defs.forEach(function (def) {
    var group = document.createElement('div');
    group.className = 'field-group';
    var label = document.createElement('label');
    label.className = 'field-label';
    label.textContent = def.label;
    label.setAttribute('for', 'param-' + def.id);
    var input = document.createElement('input');
    input.type = def.optional ? 'text' : 'number';
    input.className = 'field-input';
    input.id = 'param-' + def.id;
    input.step = def.step;
    input.value = pvals[def.id] === undefined ? '' : pvals[def.id];
    input.placeholder = def.optional ? 'e.g. 1.5' : '';
    input.addEventListener('change', function () {
      pvals[def.id] = def.optional ? input.value : parseFloat(input.value);
    });
    group.appendChild(label);
    group.appendChild(input);
    dom.paramRow.appendChild(group);
  });
}

function renderPointsTable() {
  dom.pointsTable.innerHTML = '';
  state.points.forEach(function (pt, idx) {
    var row = document.createElement('div');
    row.className = 'points-row';

    var xInput = document.createElement('input');
    xInput.type = 'number'; xInput.step = 'any'; xInput.value = pt.x;
    xInput.setAttribute('aria-label', 'x value for point ' + (idx + 1));
    xInput.addEventListener('change', function () { state.points[idx].x = parseFloat(xInput.value); });

    var yInput = document.createElement('input');
    yInput.type = 'number'; yInput.step = 'any'; yInput.value = pt.y;
    yInput.setAttribute('aria-label', 'y value for point ' + (idx + 1));
    yInput.addEventListener('change', function () { state.points[idx].y = parseFloat(yInput.value); });

    var removeBtn = document.createElement('button');
    removeBtn.type = 'button';
    removeBtn.className = 'points-row__remove';
    removeBtn.innerHTML = '&times;';
    removeBtn.setAttribute('aria-label', 'Remove point ' + (idx + 1));
    removeBtn.addEventListener('click', function () {
      if (state.points.length <= 3) return;
      state.points.splice(idx, 1);
      renderPointsTable();
    });

    row.appendChild(xInput);
    row.appendChild(yInput);
    row.appendChild(removeBtn);
    dom.pointsTable.appendChild(row);
  });
}

function switchMethod(method) {
  state.method = method;
  state.steps = [];
  state.computeError = null;
  state.currentIndex = 0;
  renderMethodTabs();
  dom.equationGroup.hidden = method === 'spline';
  dom.splinePointsGroup.hidden = method !== 'spline';
  dom.configNote.textContent = method === 'newton'
    ? 'Uses the symbolic derivative of f(x), differentiated automatically.'
    : method === 'secant'
      ? 'Approximates the derivative using the slope between two consecutive points.'
      : method === 'muller'
        ? 'Fits a parabola through three points; can locate complex roots.'
        : 'Builds a smooth, natural cubic spline through your data points.';
  renderParamRow();
  if (method === 'spline') renderPointsTable();
  validateEquation();
  renderAll();
}

function validateEquation() {
  if (state.method === 'spline') return true;
  var exprStr = dom.equationInput.value.trim();
  var parsed = compileFunction(exprStr);
  if (parsed.error) {
    dom.equationStatus.textContent = 'Parse error: ' + parsed.error;
    dom.equationStatus.classList.add('is-error');
    return false;
  }
  if (state.method === 'newton' && !parsed.fPrime) {
    dom.equationStatus.textContent = 'Parsed, but no symbolic derivative could be found.';
    dom.equationStatus.classList.add('is-error');
    return false;
  }
  dom.equationStatus.textContent = 'Parsed successfully.';
  dom.equationStatus.classList.remove('is-error');
  return true;
}

/*COMPUTE*/
function compute() {
  var result;
  if (state.method === 'spline') {
    var pts = state.points.slice().sort(function (a, b) { return a.x - b.x; });
    var xs = pts.map(function (p) { return p.x; });
    var seen = {};
    var dup = xs.some(function (v) { if (seen[v]) return true; seen[v] = true; return false; });
    if (pts.length < 3) {
      result = { steps: [], error: 'Add at least 3 data points to build a spline.' };
    } else if (dup) {
      result = { steps: [], error: 'Data points must have distinct x values.' };
    } else {
      state.points = pts;
      renderPointsTable();
      result = runSpline(pts);
    }
  } else {
    if (!validateEquation()) {
      result = { steps: [], error: 'Fix the expression above before computing.' };
    } else {
      var exprStr = dom.equationInput.value.trim();
      var p = state.params[state.method];
      if (state.method === 'newton') {
        result = runNewton(exprStr, p.x0, p.tol, Math.max(1, Math.round(p.maxIter)));
      } else if (state.method === 'secant') {
        result = runSecant(exprStr, p.x0, p.x1, p.tol, Math.max(1, Math.round(p.maxIter)));
      } else {
        result = runMuller(exprStr, p.x0, p.x1, p.x2, p.tol, Math.max(1, Math.round(p.maxIter)));
      }
    }
  }
  state.steps = result.steps;
  state.computeError = result.error;
  state.currentIndex = 0;
  stopPlaying();
  renderAll();
  triggerStepAnimation();
}

/*TRANSPORT*/
function clampIndex(i) {
  return Math.max(0, Math.min(state.steps.length - 1, i));
}
function goToStep(i) {
  var next = clampIndex(i);
  if (next === state.currentIndex) return;
  state.currentIndex = next;
  renderAll();
  triggerStepAnimation();
}
function stepNext() { goToStep(state.currentIndex + 1); }
function stepPrev() { goToStep(state.currentIndex - 1); }
function resetSteps() { goToStep(0); }

function stopPlaying() {
  state.playing = false;
  if (state.playTimer) { clearInterval(state.playTimer); state.playTimer = null; }
  dom.playIcon.innerHTML = '<path d="M8 5 L19 12 L8 19 Z" fill="currentColor"/>';
}
function togglePlay() {
  if (!state.steps.length) return;
  if (state.playing) { stopPlaying(); return; }
  state.playing = true;
  dom.playIcon.innerHTML = '<path d="M7 5 H10 V19 H7 Z M14 5 H17 V19 H14 Z" fill="currentColor"/>';
  state.playTimer = setInterval(function () {
    if (state.currentIndex >= state.steps.length - 1) { stopPlaying(); return; }
    stepNext();
  }, 900);
}

/*CANVAS SETUP*/
function setupCanvas() {
  var canvas = dom.canvas;
  var frame = canvas.parentElement;
  var rect = frame.getBoundingClientRect();
  var dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.round(rect.width * dpr));
  canvas.height = Math.max(1, Math.round(rect.height * dpr));
  canvas.style.width = rect.width + 'px';
  canvas.style.height = rect.height + 'px';
  var ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx: ctx, width: rect.width, height: rect.height };
}

function niceStep(range, targetTicks) {
  targetTicks = targetTicks || 8;
  if (!isFinite(range) || range <= 0) return 1;
  var rough = range / targetTicks;
  var mag = Math.pow(10, Math.floor(Math.log10(rough)));
  var norm = rough / mag;
  var step;
  if (norm < 1.5) step = 1; else if (norm < 3) step = 2; else if (norm < 7) step = 5; else step = 10;
  return step * mag;
}

function makeTransform(bounds, w, h) {
  var margin = { l: 54, r: 20, t: 20, b: 36 };
  var pw = Math.max(1, w - margin.l - margin.r);
  var ph = Math.max(1, h - margin.t - margin.b);
  var sx = pw / (bounds.xMax - bounds.xMin);
  var sy = ph / (bounds.yMax - bounds.yMin);
  return {
    margin: margin, pw: pw, ph: ph,
    toPx: function (x, y) {
      return [margin.l + (x - bounds.xMin) * sx, margin.t + ph - (y - bounds.yMin) * sy];
    }
  };
}

function drawGrid(ctx, bounds, t, w, h) {
  ctx.save();
  ctx.strokeStyle = 'rgba(0,0,0,0.08)';
  ctx.lineWidth = 1;
  ctx.font = '10px "Courier New", monospace';
  ctx.fillStyle = 'rgba(0,0,0,0.45)';
  var xStep = niceStep(bounds.xMax - bounds.xMin);
  var yStep = niceStep(bounds.yMax - bounds.yMin);
  var gx;
  for (gx = Math.ceil(bounds.xMin / xStep) * xStep; gx <= bounds.xMax; gx += xStep) {
    var px = t.toPx(gx, 0)[0];
    ctx.beginPath(); ctx.moveTo(px, t.margin.t); ctx.lineTo(px, h - t.margin.b); ctx.stroke();
    ctx.fillText(fmt(gx, 3), px + 3, h - t.margin.b + 15);
  }
  var gy;
  for (gy = Math.ceil(bounds.yMin / yStep) * yStep; gy <= bounds.yMax; gy += yStep) {
    var py = t.toPx(0, gy)[1];
    ctx.beginPath(); ctx.moveTo(t.margin.l, py); ctx.lineTo(w - t.margin.r, py); ctx.stroke();
    ctx.fillText(fmt(gy, 3), 6, py - 3);
  }
  ctx.strokeStyle = 'rgba(0,0,0,0.35)';
  ctx.lineWidth = 1.3;
  if (bounds.yMin <= 0 && bounds.yMax >= 0) {
    var py0 = t.toPx(0, 0)[1];
    ctx.beginPath(); ctx.moveTo(t.margin.l, py0); ctx.lineTo(w - t.margin.r, py0); ctx.stroke();
  }
  if (bounds.xMin <= 0 && bounds.xMax >= 0) {
    var px0 = t.toPx(0, 0)[0];
    ctx.beginPath(); ctx.moveTo(px0, t.margin.t); ctx.lineTo(px0, h - t.margin.b); ctx.stroke();
  }
  ctx.restore();
}

function drawCurveFromSamples(ctx, t, bounds, sampleFn, color, width) {
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = width || 2;
  ctx.beginPath();
  var N = 420, started = false;
  var span = bounds.yMax - bounds.yMin;
  for (var i = 0; i <= N; i++) {
    var x = bounds.xMin + (bounds.xMax - bounds.xMin) * i / N;
    var y = sampleFn(x);
    if (y === null || y === undefined || !isFinite(y) || y < bounds.yMin - span || y > bounds.yMax + span) {
      started = false;
      continue;
    }
    var p = t.toPx(x, y);
    if (!started) { ctx.moveTo(p[0], p[1]); started = true; } else { ctx.lineTo(p[0], p[1]); }
  }
  ctx.stroke();
  ctx.restore();
}

function drawPointMarker(ctx, t, x, y, color, r) {
  var p = t.toPx(x, y);
  ctx.beginPath();
  ctx.arc(p[0], p[1], r || 4, 0, Math.PI * 2);
  ctx.fillStyle = color;
  ctx.fill();
  ctx.lineWidth = 1.4;
  ctx.strokeStyle = 'rgba(255,255,255,0.9)';
  ctx.stroke();
}

function drawDashedSeg(ctx, t, x1, y1, x2, y2, color, alpha) {
  var p1 = t.toPx(x1, y1), p2 = t.toPx(x2, y2);
  ctx.save();
  ctx.globalAlpha = alpha === undefined ? 1 : alpha;
  ctx.strokeStyle = color; ctx.lineWidth = 1.2; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(p1[0], p1[1]); ctx.lineTo(p2[0], p2[1]); ctx.stroke();
  ctx.restore();
}

function drawSolidSeg(ctx, t, x1, y1, x2, y2, color, width, glow, alpha) {
  var p1 = t.toPx(x1, y1), p2 = t.toPx(x2, y2);
  ctx.save();
  ctx.globalAlpha = alpha === undefined ? 1 : alpha;
  if (glow) { ctx.shadowColor = color; ctx.shadowBlur = 4; }
  ctx.strokeStyle = color; ctx.lineWidth = width || 2;
  ctx.beginPath(); ctx.moveTo(p1[0], p1[1]); ctx.lineTo(p2[0], p2[1]); ctx.stroke();
  ctx.restore();
}

var COLOR_CURVE = '#3366cc';
var COLOR_PHOSPHOR = '#e08a2b';
var COLOR_DIM = '#aaaaaa';
var COLOR_DANGER = '#cc3333';

/*BOUNDS PER METHOD*/
function rootFindingBounds(f, steps) {
  var xs = [];
  steps.forEach(function (s) {
    var rx = reOf(s.x);
    if (rx !== null && isFinite(rx)) xs.push(rx);
    if (typeof s.xNext === 'number' && isFinite(s.xNext)) xs.push(s.xNext);
    if (s.quadPoints) s.quadPoints.forEach(function (q) { if (isFinite(q.re)) xs.push(q.re); });
  });
  if (!xs.length) xs = [-5, 5];
  var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
  if (xMin === xMax) { xMin -= 2; xMax += 2; }
  var pad = (xMax - xMin) * 0.35 || 2;
  xMin -= pad; xMax += pad;

  var ys = [];
  var N = 240;
  for (var i = 0; i <= N; i++) {
    var x = xMin + (xMax - xMin) * i / N;
    var y = safeEval(f, x);
    if (y !== null && isFinite(y)) ys.push(y);
  }
  if (!ys.length) ys = [-5, 5];
  ys.sort(function (a, b) { return a - b; });
  var lo = ys[Math.floor(ys.length * 0.03)];
  var hi = ys[Math.min(ys.length - 1, Math.ceil(ys.length * 0.97))];
  var yMin = Math.min(lo, 0), yMax = Math.max(hi, 0);
  if (yMin === yMax) { yMin -= 2; yMax += 2; }
  var ypad = (yMax - yMin) * 0.2 || 2;
  return { xMin: xMin, xMax: xMax, yMin: yMin - ypad, yMax: yMax + ypad };
}

function splineBounds(points, steps) {
  var xs = points.map(function (p) { return p.x; });
  var ys = points.map(function (p) { return p.y; });
  var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
  var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
  var xpad = (xMax - xMin) * 0.15 || 1;
  var ypad = (yMax - yMin) * 0.3 || 1;
  return { xMin: xMin - xpad, xMax: xMax + xpad, yMin: yMin - ypad, yMax: yMax + ypad };
}

/*METHOD OVERLAYS*/
function tangentLineY(x, x0, y0, slope) { return y0 + slope * (x - x0); }

function drawNewtonOverlay(ctx, t, bounds, steps, idx, animT) {
  for (var i = 0; i <= idx; i++) {
    var s = steps[i];
    if (s.fx === null || s.fx === undefined) continue;
    var isCurrent = i === idx;
    var alpha = isCurrent ? 1 : 0.4;
    drawDashedSeg(ctx, t, s.x, 0, s.x, s.fx, COLOR_DIM, alpha);
    drawPointMarker(ctx, t, s.x, s.fx, isCurrent ? COLOR_PHOSPHOR : '#666666', isCurrent ? 5 : 3.4);
    if (s.fpx !== undefined && s.fpx !== null && s.xNext !== undefined) {
      var ends = { y1: tangentLineY(bounds.xMin, s.x, s.fx, s.fpx), y2: tangentLineY(bounds.xMax, s.x, s.fx, s.fpx) };
      var a = isCurrent ? Math.max(0.25, animT) : alpha;
      drawSolidSeg(ctx, t, bounds.xMin, ends.y1, bounds.xMax, ends.y2, isCurrent ? COLOR_PHOSPHOR : COLOR_DIM, isCurrent ? 2.2 : 1.2, isCurrent, a);
      var revealT = isCurrent ? animT : 1;
      var mx = s.x + (s.xNext - s.x) * revealT;
      var my = s.fx + (0 - s.fx) * revealT;
      if (isCurrent) drawPointMarker(ctx, t, mx, my, COLOR_CURVE, 4.2);
      if (!isCurrent || animT > 0.85) drawPointMarker(ctx, t, s.xNext, 0, '#555555', 3);
    }
  }
}

function drawSecantOverlay(ctx, t, bounds, steps, idx, animT) {
  for (var i = 0; i <= idx; i++) {
    var s = steps[i];
    if (s.fx === null || s.fx === undefined) continue;
    var isCurrent = i === idx;
    drawPointMarker(ctx, t, s.x, s.fx, isCurrent ? COLOR_PHOSPHOR : '#666666', isCurrent ? 5 : 3.4);
  }
  for (i = 2; i <= idx; i++) {
    var p0 = steps[i - 2], p1 = steps[i - 1], target = steps[i];
    if (p0.fx === null || p1.fx === null) continue;
    var isCurrent = i === idx;
    if (p1.x === p0.x) continue;
    var slope = (p1.fx - p0.fx) / (p1.x - p0.x);
    var y1 = tangentLineY(bounds.xMin, p1.x, p1.fx, slope);
    var y2 = tangentLineY(bounds.xMax, p1.x, p1.fx, slope);
    var a = isCurrent ? Math.max(0.25, animT) : 0.4;
    drawSolidSeg(ctx, t, bounds.xMin, y1, bounds.xMax, y2, isCurrent ? COLOR_PHOSPHOR : COLOR_DIM, isCurrent ? 2.2 : 1.2, isCurrent, a);
    if (target && target.x !== undefined && isFinite(target.x)) {
      var revealT = isCurrent ? animT : 1;
      var mx = p1.x + (target.x - p1.x) * revealT;
      var my = p1.fx + (0 - p1.fx) * revealT;
      if (isCurrent) drawPointMarker(ctx, t, mx, my, COLOR_CURVE, 4.2);
    }
  }
}

function lagrangeQuad(p0, p1, p2) {
  return function (x) {
    var l0 = ((x - p1.x) * (x - p2.x)) / ((p0.x - p1.x) * (p0.x - p2.x));
    var l1 = ((x - p0.x) * (x - p2.x)) / ((p1.x - p0.x) * (p1.x - p2.x));
    var l2 = ((x - p0.x) * (x - p1.x)) / ((p2.x - p0.x) * (p2.x - p1.x));
    return p0.y * l0 + p1.y * l1 + p2.y * l2;
  };
}

function drawMullerOverlay(ctx, t, bounds, steps, idx, animT) {
  for (var i = 0; i <= idx; i++) {
    var s = steps[i];
    if (!s.quadPoints) continue;
    var isCurrent = i === idx;
    var qp = s.quadPoints, qf = s.quadFs;
    var allReal = qp.every(function (p) { return Math.abs(p.im) < 1e-6; }) && qf.every(function (p) { return Math.abs(p.im) < 1e-6; });
    var alpha = isCurrent ? 1 : 0.35;
    if (allReal) {
      var real0 = { x: qp[0].re, y: qf[0].re };
      var real1 = { x: qp[1].re, y: qf[1].re };
      var real2 = { x: qp[2].re, y: qf[2].re };
      var quadFn = lagrangeQuad(real0, real1, real2);
      ctx.save();
      ctx.globalAlpha = alpha;
      drawCurveFromSamples(ctx, t, bounds, quadFn, isCurrent ? COLOR_PHOSPHOR : COLOR_DIM, isCurrent ? 2.1 : 1.2);
      ctx.restore();
      [real0, real1, real2].forEach(function (p) {
        drawPointMarker(ctx, t, p.x, p.y, isCurrent ? '#555555' : '#999999', isCurrent ? 3.6 : 2.8);
      });
      if (Math.abs(s.x.im) < 1e-6) {
        var revealT = isCurrent ? animT : 1;
        var startY = quadFn(real2.x);
        var mx = real2.x + (s.x.re - real2.x) * revealT;
        var my = startY + (0 - startY) * revealT;
        if (isCurrent) drawPointMarker(ctx, t, mx, my, COLOR_CURVE, 4.2);
        if (!isCurrent || animT > 0.85) drawPointMarker(ctx, t, s.x.re, 0, '#555555', 3);
      }
    }
  }
}

function drawSplineOverlay(ctx, t, bounds, steps, idx, animT, points) {
  points.forEach(function (p) {
    drawPointMarker(ctx, t, p.x, p.y, '#777777', 3.6);
  });
  for (var i = 0; i <= idx; i++) {
    var seg = steps[i];
    var isCurrent = i === idx;
    var reveal = isCurrent ? Math.max(0.02, animT) : 1;
    var xEnd = seg.x0 + (seg.x1 - seg.x0) * reveal;
    var fn = function (x) {
      var dx = x - seg.x0;
      return seg.a + seg.b * dx + seg.c * dx * dx + seg.d * dx * dx * dx;
    };
    ctx.save();
    var clipBounds = { xMin: seg.x0, xMax: xEnd, yMin: bounds.yMin, yMax: bounds.yMax };
    if (xEnd > seg.x0) drawCurveFromSamples(ctx, t, clipBounds, fn, isCurrent ? COLOR_CURVE : COLOR_DIM, isCurrent ? 2.6 : 1.6);
    ctx.restore();
    drawPointMarker(ctx, t, seg.x0, seg.a, isCurrent ? COLOR_PHOSPHOR : '#666666', isCurrent ? 4.6 : 3.2);
    if (!isCurrent || animT > 0.9) {
      drawPointMarker(ctx, t, seg.x1, fn(seg.x1), isCurrent ? COLOR_PHOSPHOR : '#666666', isCurrent ? 4.6 : 3.2);
    }
  }
  var evalXRaw = state.params.spline.evalX;
  var evalX = parseFloat(evalXRaw);
  if (evalXRaw !== '' && isFinite(evalX) && idx === steps.length - 1) {
    var y = evalSplineAt(steps, evalX);
    if (y !== null) {
      drawDashedSeg(ctx, t, evalX, bounds.yMin, evalX, y, COLOR_CURVE, 0.8);
      drawDashedSeg(ctx, t, bounds.xMin, y, evalX, y, COLOR_CURVE, 0.8);
      drawPointMarker(ctx, t, evalX, y, COLOR_CURVE, 5.5);
    }
  }
}

/*MAIN CANVAS RENDER*/
function renderCanvas() {
  var setup = setupCanvas();
  var ctx = setup.ctx, w = setup.width, h = setup.height;
  ctx.clearRect(0, 0, w, h);
  if (!state.steps.length) return;

  var idx = state.currentIndex;
  var animT = state.animT;

  if (state.method === 'spline') {
    var bounds = splineBounds(state.points, state.steps);
    var t = makeTransform(bounds, w, h);
    drawGrid(ctx, bounds, t, w, h);
    drawSplineOverlay(ctx, t, bounds, state.steps, idx, animT, state.points);
  } else {
    var exprStr = dom.equationInput.value.trim();
    var parsed = compileFunction(exprStr);
    if (parsed.error) return;
    var bounds2 = rootFindingBounds(parsed.f, state.steps);
    var t2 = makeTransform(bounds2, w, h);
    drawGrid(ctx, bounds2, t2, w, h);
    drawCurveFromSamples(ctx, t2, bounds2, function (x) { return safeEval(parsed.f, x); }, COLOR_CURVE, 2);
    if (state.method === 'newton') drawNewtonOverlay(ctx, t2, bounds2, state.steps, idx, animT);
    else if (state.method === 'secant') drawSecantOverlay(ctx, t2, bounds2, state.steps, idx, animT);
    else if (state.method === 'muller') drawMullerOverlay(ctx, t2, bounds2, state.steps, idx, animT);

    var cur = state.steps[idx];
    if (cur && cur.failed) {
      var rx = reOf(cur.x);
      if (rx !== null && isFinite(rx)) {
        var py = cur.fx !== null && cur.fx !== undefined ? reOf(cur.fx) : 0;
        var p = t2.toPx(rx, py || 0);
        ctx.save();
        ctx.strokeStyle = COLOR_DANGER;
        ctx.lineWidth = 2;
        var pulse = 6 + 4 * Math.sin(performance.now() / 220);
        ctx.beginPath(); ctx.arc(p[0], p[1], pulse, 0, Math.PI * 2); ctx.stroke();
        ctx.restore();
      }
    }
  }
}

function prefersReducedMotion() {
  return window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
}

function triggerStepAnimation() {
  if (state.raf) cancelAnimationFrame(state.raf);
  if (prefersReducedMotion()) { state.animT = 1; renderCanvas(); return; }
  state.animT = 0;
  var start = performance.now();
  var dur = 420;
  function frame(now) {
    var elapsed = now - start;
    var lt = Math.min(1, elapsed / dur);
    state.animT = 1 - Math.pow(1 - lt, 3);
    renderCanvas();
    if (lt < 1) {
      state.raf = requestAnimationFrame(frame);
    } else if (state.steps[state.currentIndex] && state.steps[state.currentIndex].failed) {
      state.raf = requestAnimationFrame(pulseFrame);
    }
  }
  function pulseFrame() {
    renderCanvas();
    state.raf = requestAnimationFrame(pulseFrame);
  }
  state.raf = requestAnimationFrame(frame);
}

/*TABLE RENDERING*/
var tableColumns = {
  newton: [
    { key: 'n', label: 'n' },
    { key: 'x', label: 'x\u2099' },
    { key: 'fx', label: 'f(x\u2099)' },
    { key: 'fpx', label: "f'(x\u2099)" },
    { key: 'relError', label: 'Rel. error %' }
  ],
  secant: [
    { key: 'n', label: 'n' },
    { key: 'x', label: 'x\u2099' },
    { key: 'fx', label: 'f(x\u2099)' },
    { key: 'relError', label: 'Rel. error %' }
  ],
  muller: [
    { key: 'n', label: 'n' },
    { key: 'x', label: 'x\u2099' },
    { key: 'fx', label: '|f(x\u2099)|' },
    { key: 'relError', label: 'Rel. error %' }
  ],
  spline: [
    { key: 'n', label: 'segment' },
    { key: 'interval', label: 'interval [x\u2c7c, x\u2c7c\u208a\u2081]' },
    { key: 'a', label: 'a' },
    { key: 'b', label: 'b' },
    { key: 'c', label: 'c' },
    { key: 'd', label: 'd' }
  ]
};

function renderTable() {
  var cols = tableColumns[state.method];
  dom.dataTableHead.innerHTML = '<tr>' + cols.map(function (c) { return '<th>' + c.label + '</th>'; }).join('') + '</tr>';
  dom.dataTableBody.innerHTML = '';
  if (!state.steps.length) {
    dom.dataTable.classList.add('hidden');
    dom.tableEmpty.classList.remove('hidden');
    return;
  }
  dom.dataTable.classList.remove('hidden');
  dom.tableEmpty.classList.add('hidden');

  for (var i = 0; i <= state.currentIndex; i++) {
    var s = state.steps[i];
    var tr = document.createElement('tr');
    if (s.failed) tr.classList.add('is-error-row');
    var cells = cols.map(function (c) {
      if (state.method === 'spline') {
        if (c.key === 'interval') return '[' + fmt(s.x0, 4) + ', ' + fmt(s.x1, 4) + ']';
        if (c.key === 'n') return String(s.n);
        return fmt(s[c.key], 5);
      }
      if (c.key === 'n') return String(s.n);
      if (c.key === 'x') return fmtComplex(s.x);
      if (c.key === 'fx') {
        if (s.fx === null || s.fx === undefined) return '\u2014';
        if (state.method === 'muller') return fmt(math.abs(s.fx), 5);
        return fmt(s.fx, 6);
      }
      if (c.key === 'fpx') return s.fpx === undefined || s.fpx === null ? '\u2014' : fmt(s.fpx, 6);
      if (c.key === 'relError') return s.relError === undefined ? '\u2014' : fmt(s.relError, 5) + '%';
      return '\u2014';
    });
    tr.innerHTML = cells.map(function (v) { return '<td>' + v + '</td>'; }).join('');
    dom.dataTableBody.appendChild(tr);
  }
}

/*READOUT / TRANSPORT / ERROR CARD*/
function renderReadout() {
  dom.readoutMethod.textContent = methodNames[state.method];
  var cur = state.steps[state.currentIndex];
  if (!state.steps.length) {
    dom.readoutState.textContent = 'READY';
    dom.readoutState.className = 'readout__value';
    dom.readoutRoot.textContent = '\u2014';
    return;
  }
  if (cur && cur.failed) {
    dom.readoutState.textContent = 'ERROR';
    dom.readoutState.className = 'readout__value is-error';
  } else if (cur && cur.converged) {
    dom.readoutState.textContent = 'CONVERGED';
    dom.readoutState.className = 'readout__value is-success';
  } else if (state.currentIndex === state.steps.length - 1 && cur && cur.maxedOut) {
    dom.readoutState.textContent = 'MAX ITER';
    dom.readoutState.className = 'readout__value';
  } else {
    dom.readoutState.textContent = 'STEPPING';
    dom.readoutState.className = 'readout__value';
  }
  if (state.method === 'spline') {
    dom.readoutRoot.textContent = 'n/a';
  } else {
    var lastGood = null;
    for (var i = state.currentIndex; i >= 0; i--) {
      if (!state.steps[i].failed) { lastGood = state.steps[i]; break; }
    }
    dom.readoutRoot.textContent = lastGood ? fmtComplex(lastGood.x, 5) : '\u2014';
  }
}

function renderErrorCard() {
  var cur = state.steps[state.currentIndex];
  if (cur && cur.failed) {
    dom.errorCard.classList.remove('hidden');
    dom.errorCardTitle.textContent = 'Method breakdown at step ' + cur.n;
    dom.errorCardBody.textContent = cur.errorMsg || 'The iteration could not continue from this point.';
  } else {
    dom.errorCard.classList.add('hidden');
  }
}

function renderTransport() {
  var total = state.steps.length;
  dom.stepCurrent.textContent = total ? String(state.currentIndex + 1) : '0';
  dom.stepTotal.textContent = String(total);
  dom.stepSlider.max = String(Math.max(0, total - 1));
  dom.stepSlider.value = String(state.currentIndex);
  dom.stepSlider.disabled = total <= 1;

  var atStart = state.currentIndex <= 0;
  var atEnd = !total || state.currentIndex >= total - 1;
  dom.btnPrev.disabled = atStart || !total;
  dom.btnNext.disabled = atEnd;
  dom.btnReset.disabled = !total;
  dom.btnPlay.disabled = !total || atEnd && !state.playing;

  if (atEnd) stopPlaying();

  var legendItems = [];
  if (state.method === 'newton') {
    legendItems = [['tangent line', COLOR_PHOSPHOR], ['f(x) curve', COLOR_CURVE], ['next estimate', '#555555']];
  } else if (state.method === 'secant') {
    legendItems = [['secant line', COLOR_PHOSPHOR], ['f(x) curve', COLOR_CURVE]];
  } else if (state.method === 'muller') {
    legendItems = [['interpolating parabola', COLOR_PHOSPHOR], ['f(x) curve', COLOR_CURVE]];
  } else {
    legendItems = [['current segment', COLOR_CURVE], ['completed segments', COLOR_DIM], ['data points', '#777777']];
  }
  dom.transportLegend.innerHTML = legendItems.map(function (it) {
    return '<span><span class="legend-swatch" style="background:' + it[1] + '"></span>' + it[0] + '</span>';
  }).join('');

  dom.canvasEmpty.style.display = total ? 'none' : 'flex';
}

function renderConfigError() {
  if (state.computeError) {
    dom.equationStatus.textContent = state.computeError;
    dom.equationStatus.classList.add('is-error');
  }
}

function renderAll() {
  renderReadout();
  renderTransport();
  renderErrorCard();
  renderTable();
  renderConfigError();
  renderCanvas();
}

/*EVENTS*/
function wireEvents() {
  dom.methodTabs.addEventListener('click', function (e) {
    var btn = e.target.closest('.method-tab');
    if (!btn) return;
    switchMethod(btn.getAttribute('data-method'));
  });

  dom.equationInput.addEventListener('input', function () { validateEquation(); });

  dom.addPointBtn.addEventListener('click', function () {
    var lastX = state.points.length ? state.points[state.points.length - 1].x : 0;
    state.points.push({ x: lastX + 1, y: 0 });
    renderPointsTable();
  });

  dom.computeBtn.addEventListener('click', compute);

  dom.btnReset.addEventListener('click', function () { stopPlaying(); resetSteps(); });
  dom.btnPrev.addEventListener('click', function () { stopPlaying(); stepPrev(); });
  dom.btnNext.addEventListener('click', function () { stopPlaying(); stepNext(); });
  dom.btnPlay.addEventListener('click', togglePlay);
  dom.stepSlider.addEventListener('input', function () {
    stopPlaying();
    goToStep(parseInt(dom.stepSlider.value, 10));
  });

  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(renderCanvas, 80);
  });
}

/*INIT*/
function init() {
  cacheDom();
  renderMethodTabs();
  renderParamRow();
  wireEvents();
  validateEquation();
  compute();
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', init);
} else {
  init();
}
