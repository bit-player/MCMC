
// Explore variants of the Metropolis algorithm
// for the Ising model of ferromagnetism. The variants
// differ in the order in which sites are visited during
// each sweep through the lattice.

// MIT License   Copyright (c) 2021 Brian Hayes

// This program accompanies the article "Three Months in Monte Carlo"
// at http://bit-player.org/2021/three-months-in-monte-carlo


// Get references to HTML elements, link event handlers

const theViz = document.getElementById("the-viz");
const theCanvas = document.getElementById("the-canvas");
const ctx = theCanvas.getContext("2d");

const Mreadout = document.getElementById("magnetization-readout");
const Rreadout = document.getElementById("correlation-readout");

const theRunButton = document.getElementById("run-button");
theRunButton.onclick = doRunButton;
const theStepButton = document.getElementById("step-button");
theStepButton.onclick = doStepButton;
const theResetButton = document.getElementById("reset-button");
theResetButton.onclick = doResetButton;

const theTemperatureSlider = document.getElementById("temperature-slider");
theTemperatureSlider.onchange = adjustTemperature;
theTemperatureSlider.oninput = adjustTemperature;
const temperatureReadout = document.getElementById("temperature-readout");

const radioButtons = document.querySelectorAll('input[name="visit-sequence"]');
for (rb of radioButtons) {
  rb.addEventListener('change', chooseAlgorithm);
}

const outputDiv = document.getElementById("stats-output");


// Set attributes of the lattice display

const gridSize = 100;
const N = gridSize * gridSize;
const gridEdge = gridSize - 1;
const cellSize = theCanvas.width / gridSize;
const upColor = "#4b0082";           // indigo
const downColor = "#a297ff";         // call it mauve, though it's really light blue
let temperature = Number(theTemperatureSlider.value);


// Global variables for the state machine controlling program execution

let timer;
let algorithm = 'typewriter'   
let state = 'paused';
let thePermutation = randomPermutation(N);


// An association list connecting strings to functions,
// and also the speeds at which those functions run.

const algoList = 
  {'typewriter': {'fn': runMetro_typewriter, 'speed': 30},
   'checkerboard': {'fn': runMetro_checkerboard, 'speed': 30},
   'diagonal': {'fn': runMetro_diagonal, 'speed': 30},
   'boustrophedon': {'fn': runMetro_boustrophedon, 'speed': 30},
   'permuted': {'fn': runMetro_permuted, 'speed': 30},
   'repermuted': {'fn': runMetro_re_permuted, 'speed': 30},
   'random': {'fn': runMetro_random, 'speed': 30},
   'simultaneous': {'fn': runMetro_simultaneous, 'speed': 200}
  }



// build the lattice as an array of arrays of cells

let lattice = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    lattice[i] = new Array(gridSize);
}


// The 'simultaneous' visitation order needs a place to stash
// information on which spins are to be flipped, once the sweep
// of the lattice is complete.

let shadowLattice = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    shadowLattice[i] = new Array(gridSize);
}


// Update the canvas and the two meters. Called after
// every macrostep.

function drawLattice() {
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      ctx.fillStyle = (lattice[x][y] === 1) ? upColor : downColor;
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
  Mreadout.innerHTML = formatReadout(magnetization());
  Rreadout.innerHTML = formatReadout(local_correlations());
}


// The following pair of functions does the calculation of
// magnetization and local-correlation values. M is just the
// sum of the N +1 and -1 spin values, divided by N to
// yield a value between -1 and +1. 

// R is defined as the number of like neighboring pairs minus
// the number of unlike pairs. But rather than write a function
// to do that calculation, we can use the existing function
// calcDeltaE, which returns 8 times the desired sum.

function magnetization() {
  let M = 0;
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      M += lattice[x][y];
    }
  }
  return M / N;  
}

function local_correlations() {
  let R = 0;
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      R += calcDeltaE(x, y) / 8;
    }
  }
  return R / N;
}


// To avoid jitter on the screen, always show a sign prefix,
// and always show three decimal places. (Also use monospace
// font, but that's specified in the CSS.)

function formatReadout(val) {
  sign = (val < 0) ? '&minus;' : '&plus;'
  return sign + Math.abs(val).toFixed(3);
}


// Used during initialization to produce random pattern.

function coinFlip() {
  return Math.random() < 0.5;
}


// The 'permuted' and 're-permuted' visitation sequences need
// random permutations of the sequence [0..n-1]

function randomPermutation(n) {
  const p = Array.from({length: n}, (_, i) => i);
  let temp;
  for (let i = n-1; i > 1; i--) {
    let j = Math.floor(Math.random() * n);
    temp = p[j];
    p[j] = p[i];
    p[i] = temp;
  }
  return p;
}

// Create and display an initial, fully randomized array of cells

function init() {
  let x, y;
  state = 'scrambled';
  thePermutation = randomPermutation(N);
  for (y = 0; y < gridSize; y++) {
    for (x = 0; x < gridSize; x++) {
      lattice[x][y] = coinFlip() ? 1 : -1;
    }
  }
  drawLattice();
}


// Given the coordinates of a cell, survey its four neighbors
// and compute the sum of their spins, which will be a value
// in the set {-4, -2, 0, 2, 4}. Multiply by twice the value of the
// central spin, yielding an answer in {-8, -4, 0, 4, 8}. This
// is the net change in system energy if the spin were to flip.

function calcDeltaE(x, y) {
  let north, south, east, west, deltaE, boltz;
  north = lattice[x][(y > 0) ? y - 1 : gridEdge];   // using toroidal lattice
  south = lattice[x][(y < gridEdge) ? y + 1 : 0];
  east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
  west  = lattice[(x < gridEdge) ? x + 1 : 0][y];
  deltaE = 2 * lattice[x][y] * (north + south + east + west);
  return deltaE;
}


// Implements the M_rule acceptance function.

function updateSpin(x, y) {
  let north, south, east, west, deltaE, boltzmann;
  north = lattice[x][(y > 0) ? y - 1 : gridEdge];
  south = lattice[x][(y < gridEdge) ? y + 1 : 0];
  east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
  west  = lattice[(x < gridEdge) ? x + 1 : 0][y];
  deltaE = 2 * lattice[x][y] * (north + south + east + west);
  boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
  if ((deltaE <= 0) || (Math.random() < boltzmann)) {
    lattice[x][y] *= -1;
  }
}


// The following series of functions implements the various visitation
// sequences. Each one runs a sequence of N visits to lattice sites,
// performs the apprpriate spin updates, and then redraws the lattice.

// Typewriter: left to right and top to bottom.

function runMetro_typewriter() {
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      updateSpin(x, y);
    }
  }
  drawLattice();
}

// Checkerboard: visit all the black squares of the
// checkboard, then all the white squares.

function runMetro_checkerboard() {
  for (let x = 0; x < gridSize; x++) {
    for (let y = (x % 2); y < gridSize; y += 2) {
      updateSpin(x, y);
    }
  }
  for (let x = 0; x < gridSize; x++) {
    for (let y = ((x + 1) % 2); y < gridSize; y += 2) {
      updateSpin(x, y);
    }
  }
  drawLattice();
}

// Diagonal: A Cantor-like sequence, starting at the
// northwest corner and then following successive 
// southwest-to-northeast diagonals. First comes a 
// generator function, which produces the sequence
// of x,y coordinates, and then a consumer of that
// sequence. Note the asterisk: function*, the marker
// of a generator.

function* diagonalCoords() {
  let x = 0, y = 0, count = 0;
  while (count < N) {
    yield {x, y};
    count++
    if (x === gridEdge) {
      x = y + 1;
      y = gridEdge;
      }
    else if (y === 0) {
      y = x + 1;
      x = 0;
    }
    else {
      x += 1;
      y -= 1;
    }
  }
}

function runMetro_diagonal() {
  let diagonalIterator = diagonalCoords();
  for (xy of diagonalIterator) {
    updateSpin(xy.x, xy.y);
  }
  drawLattice();
}


// A test function that makes it a little clearer what's going
// on with the generator. 

function testDiagonal() {
  let iter = diagonalCoords();  // 'iter' is bound to an object with a 'next()'
  let xy = iter.next();         //    method. Calling 'next()' returns an object
  while (!xy.done) {            //    with two fields, 'value' and 'done'. 
    console.log(xy.value);
    xy = iter.next();
  }
}


// The boustrophedonic sequence: left to right, then right
// to left, as the ox plows a field. Again a generator
// produces the sequence of coordinates, which we iterate
// through in runMetro_boustrophedon.

function* boustrophedonCoords() {
  
  function evenP(n) {
    return (n % 2) === 0;
  }
  
  let x = 0, y = 0;
  while (y < gridSize) {
    yield {x, y};
    if (evenP(y)) {
      if (x < gridEdge) {
        x++;
      }
      else {
        y++
      }
    }
    else {
      if (x > 0) {
        x--
      }
      else {
        y++
      }
    }
  }
}

function runMetro_boustrophedon() {
  let bousIterator = boustrophedonCoords();
  for (xy of bousIterator) {
    updateSpin(xy.x, xy.y);
  }
  drawLattice();
}



// Simultaneous order can't use the updateSpin(x, y) function, because
// we don't want to do any updates until the whole show is over. So
// instead we make all flips in the shadowlattice, then make a separate
// sweep at the end, copying the shadow values back to the displayed
// lattice

function runMetro_simultaneous() {
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      shadowLattice[x][y] = lattice[x][y];
      let deltaE = calcDeltaE(x, y);
      let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
      if ((deltaE <= 0) || (Math.random() < boltzmann)) {
        shadowLattice[x][y] *= -1;
      }
    }
  }
  for (x = 0; x < gridSize; x++) {
    for (y = 0; y < gridSize; y++) {
      lattice[x][y] = shadowLattice[x][y];
    }
  }
  drawLattice();
}


// Generate N x, y pairs.

function runMetro_random() {
  for (let i = 0; i < N; i++) {
    let x = Math.floor(Math.random() * gridSize);
    let y = Math.floor(Math.random() * gridSize);
    updateSpin(x, y);
  }
  drawLattice();
}


// The permutation is random, but this one differs from 
// runMetro_random in that every site is visited exactly once,
// just in randomized order. The permutation is fixed from one
// macrostep to the next. A new permutation is generated only
// on Reset (or reloading the program).

// 'thePermutation' is an array of N integers, thus for a 
// 100 x 100 grid it covers the range from 0 to 9999. We
// break each such number down into x and y coordinates
// in [0, 99] by divide and remainder operations.

function runMetro_permuted() {
  for (let r of thePermutation) {
    let x = Math.floor(r / gridSize);
    let y = Math.floor(r % gridSize);
    updateSpin(x, y);
  }
  drawLattice();
}


// Identical to 'permuted', except we generate a new permutation
// before each macrostep.

function runMetro_re_permuted() {
  const p = randomPermutation(N);
  for (let r of p) {
    let x = Math.floor(r / gridSize);
    let y = Math.floor(r % gridSize);
    updateSpin(x, y);
  }
  drawLattice();
}



// Event handler for presses of the Run/Stop button. This
// is a toggle switch: It starts the simulation if the
// program is idle, but if the simulation is already running,
// the button stops it.

// The visitation sequence is chosen by indirection through the
// algolist array, which stores objects with the field 'fn'
// and 'speed'. 

function doRunButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    this.innerHTML = "Run";
  }
  else {
    state = 'running';
    this.innerHTML = "Stop";
    alg = algoList[algorithm];
    timer = setInterval(alg.fn, alg.speed);
  }
}


// Event handler for the single-step button. Runs the
// selected algorithm for one sweep of the lattice.

function doStepButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    theRunButton.innerHTML = "Run";
  }
  else {
    alg = algoList[algorithm];
    alg.fn();
  }
}


// Event handler for the Reset button. Clears the timer
// to stop the animation, then re-inits the program.
// Also clears the outputDiv, about which more below,
// under console-only functions. And creates a new
// random permutation for use by the 'permute' sequence.

function doResetButton(e) {
  state = 'paused';
  theRunButton.innerHTML = "Run";
  clearInterval(timer);
  thePermutation = randomPermutation(N);
  outputDiv.textContent = "";
  init();
}


// Event handler for the temperature slider. Update the
// 'temperature' global, and the readout for the slider.

function adjustTemperature(e) {
  temperature = Number(this.value);
  temperatureReadout.textContent = temperature.toFixed(2);
}



// Event handler for the group of radio buttons that
// select visitation sequence. If the
// simulation is already running, we put the change into
// effect by clearing the current timer and calling
// doRunButton() to start a new one.

function chooseAlgorithm(e) {
  algorithm = this.value;
  if (state === 'running') {
    state = 'switchingAlgorithms';
    clearInterval(timer);
    doRunButton();
  }
}
  

// During program loading, at this point we call the init()
// function. Once it returns we're set to go.

init();

// ******* END OF INTERACTIVE SIMULATION CODE *******
// ******* All functions below this line must *******
// ******* be invoked from the console.       *******






// Convert a number to a string, replacing
// any decimal point with "p", to make it suitable
// as a component of a file name.
function pointless(num) {
  return String(num).replace(".", "p");
}

// Record the progress of accommodation to an abrupt change in
// temperature. Generated data for Figure 15.

function stepResponse(algo, T1, T2, T1_steps, T2_steps, burnin, reps) {
  const prefix = `${algo}_T${pointless(T1)}_T${pointless(T2)}_steps${T1_steps + T2_steps}_reps${reps}`
  let cachedDrawLattice = drawLattice;
  drawLattice = function(i, j) { ; };
  let rs = new Array(T1_steps + T2_steps).fill(0);
  let Hs = new Array(T1_steps + T2_steps).fill(0);
  let fn = algoList[algo].fn;
  init();
  for (let i = 0; i < reps; i++) {
    temperature = T1;
    for (let w = 0; w < burnin; w++) {
      fn();
    }
    for (let j = 0; j < T1_steps; j++) {
      fn();
      rs[j] += local_correlations();
      Hs[j] += Math.abs(magnetization());
    }
    temperature = T2;
    for (let k = 0; k < T2_steps; k++) {
      fn();
      rs[T1_steps + k] += local_correlations();
      Hs[T1_steps + k] += Math.abs(magnetization());
    }
  }
  drawLattice = cachedDrawLattice;
  let r_results = rs.map(x => (x / reps).toFixed(4));
  let H_results = Hs.map(x => (x / reps).toFixed(4));
  outputDiv.textContent = prefix + "_r = [";
  for (let i = 0; i < r_results.length - 1; i++) {
    outputDiv.textContent += r_results[i] + ", ";
  }
  outputDiv.textContent += r_results[r_results.length - 1] + "]";
  outputDiv.textContent += "\n\n";
  outputDiv.textContent += prefix + "_H = [";
  for (let i = 0; i < H_results.length - 1; i++) {
    outputDiv.textContent += H_results[i] + ", ";
  }
  outputDiv.textContent += H_results[H_results.length - 1] + "]";
}

