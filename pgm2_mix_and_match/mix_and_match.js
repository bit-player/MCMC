
// The Metropolis algorithm and Glauber dynamics are
// distinguished by two main features: The order in which
// sites in the lattice are visited and the rule used to
// decide whether or not the spin at a visited site should
// be flipped. These programs separate the visitation
// sequence and acceptance function components, allowing
// you to explore all four possible combinations.

// MIT License   Copyright (c) 2021 Brian Hayes

// This program accompanies the article "Three Months in Monte Carlo"
// at http://bit-player.org/2021/three-months-in-monte-carlo


// Get references to HTML elements; link event handlers

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

const visitationButtons = document.querySelectorAll('input[name="visitation"]');
for (vb of visitationButtons) {
  vb.addEventListener('change', chooseVisitationFn);
}

const acceptanceButtons = document.querySelectorAll('input[name="acceptance-function"]');
for (ab of acceptanceButtons) {
  ab.addEventListener('change', chooseAcceptanceFn);
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
let visitationFn = 'deterministic';  // alternative is 'random'
let acceptanceFn = 'M_rule';    // alternatives are 'G_rule' and 'M_star_rule'
let state = 'paused';
let statsFlag = false;


// build the lattice as an array of arrays of cells

let lattice = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    lattice[i] = new Array(gridSize);
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

function local_correlations() {
  let R = 0;
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      R += calcDeltaE(x, y) / 8;
    }
  }
  return R / N;
}

function magnetization() {
  let M = 0;
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      M += lattice[x][y];
    }
  }
  return M / N;  
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



// Create and display an initial array of cells. Only
// the 'scrambled' option is available via the point-and-click
// interface, but the others can be invoked from the console.

function init(pattern='scrambled') {
  let i, j;
  state = 'scrambled';  
  switch (pattern) {
    case 'scrambled': {
      for (i = 0; i < gridSize; i++) {
        for (j = 0; j < gridSize; j++) {
          lattice[i][j] = coinFlip() ? 1 : -1;
        }
      }
    }
    break;
    case 'plus': {
      for (i = 0; i < gridSize; i++) {
        for (j = 0; j < gridSize; j++) {
          lattice[i][j] = 1;
        }
      }
    }
    break;
    case 'minus': {
      for (i = 0; i < gridSize; i++) {
        for (j = 0; j < gridSize; j++) {
          lattice[i][j] = -1;
        }
      }
    }
    break;
    default: console.log(`Unknown initState ${pattern}`);
  }
  drawLattice();
}



// Given the coordinates of a cell, survey its four neighbors
// and compute the sum of their spins, which will be a value
// in the set {-4, -2, 0, 2, 4}. Multiply by twice the value of the
// central spin, yielding an answer in {-8, -4, 0, 4, 8}. This
// is the net change in system energy if the spin were to flip.

// Wraparound (or toroidal or periodic) boundary conditions
// are hard-wired here. See Program 4 for another approach.

function calcDeltaE(x, y) {
  let north, south, east, west, deltaE;
  north = lattice[x][(y > 0) ? y - 1 : gridEdge];   // using toroidal lattice
  south = lattice[x][(y < gridEdge) ? y + 1 : 0];
  east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
  west  = lattice[(x < gridEdge) ? x + 1 : 0][y];
  deltaE = 2 * lattice[x][y] * (north + south + east + west);
  return deltaE;
}


// The architecture of these basic procedures for running the simulation
// is a little different here from what it was in Program 1. It will
// change again in subsequent programs.

// Run one macrostep (consisting of N microsteps), using either a
// deterministic or a random visitation order. Call 'updateSpin'
// for each site visited. At the end of the sweep, redraw the lattice.

function runBatch() {
  if (visitationFn === "deterministic") {
    for (let y = 0; y < gridSize; y++) {
      for (let x = 0; x < gridSize; x++) {
        updateSpin(x, y);
      }
    }
  }
  else if (visitationFn === "random") {
    for (let i = 0; i < N; i++) {
      let x = Math.floor(Math.random() * gridSize);
      let y = Math.floor(Math.random() * gridSize);
      updateSpin(x, y);
    }
  }
  drawLattice();
}



// Implementing the acceptance function -- all three in one hunk
// of code. First calculate deltaE, the amount by which the system's
// energy would change if the selected spin were to flip. This will
// be a number in the set {-8, -4, 0, 4, 8}. Then calculate the
// Boltzmann weight, which depends on deltaE and the temperature.
// Given these values, apply one of the three acceptance rules
// and either flip the spin (by negating its value) or not.

function updateSpin(x, y) {
  
  let deltaE = calcDeltaE(x, y);
  let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
  
  // The M_rule is the usual one for the Metropolis algorithm. It is
  // a piecewise function: If deltaE is negative or zero, we always
  // flip the spin. If deltaE is positive, we flip if a random variate
  // is less than the Boltzmann weight.
  
  // Note: Everytime I take a fresh look at this implementation, I think
  // I've made a terrible mistake. Math.random() returns a float in the
  // interval [0, 1), but the Boltzmann weight can range from 0 to positive
  // infinity. Thus it seems the proposed flip will almost always be
  // accepted. This is wrong. The else if branch of the M_rule is taken
  // only if deltaE is greater than 0, which means that -deltaE / T
  // must be negative, and so the Boltzmann weight exp(-deltaE / T)
  // will be less than 1.
  
  if (acceptanceFn === "M_rule") {
    if (deltaE <= 0) {
      lattice[x][y] *= -1;
    }
    else if (Math.random() < boltzmann) {   // deltaE > 0
        lattice[x][y] *= -1;
    }
  }
  
  // The G_rule (typically part of Glauber dynamics) flips
  // the spin is the random float lies below the curve
  // boltzmann / (1 + boltzmann).
  
  else if (acceptanceFn === "G_rule") {
    if (Math.random() < (boltzmann / (1 + boltzmann))) {
      lattice[x][y] *= -1;
    }
  }
  
  // The M_star_rule was invented to satisfy my own curiosity.
  // The oddest feature of the M_rule is *always* flipping a 
  // neutral cell. In M_star we make that case a coin flip.
  // (Turns out not to be all that interesting.)
  
  else if (acceptanceFn === "M_star_rule") {
    if (deltaE < 0) {
      lattice[x][y] *= -1;
    }
    else if (deltaE === 0) {
      if (coinFlip()) {
        lattice[x][y] *= -1;
      }
    } 
    else if (Math.random() < boltzmann) {   // deltaE > 0
        lattice[x][y] *= -1;
    }
  }
}
  
  
// Event handler for the Run/Stop button. Set a timer
// to run a batch every 30 milliseconds.

function doRunButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    this.innerHTML = "Run";
  }
  else {
    state = 'running';
    this.innerHTML = "Stop";
    timer = setInterval(runBatch, 30);
  }
}

// And for the single-step button.

function doStepButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    theRunButton.innerHTML = "Run";
  }
  else {
    runBatch();
  }
}


// And the Reset button.

function doResetButton(e) {
  state = 'paused';
  theRunButton.innerHTML = "Run";
  clearInterval(timer);
  outputDiv.textContent = "";
  init();
}


// And the slide that sets the temperature.

function adjustTemperature(e) {
  temperature = Number(this.value);
  temperatureReadout.textContent = temperature.toFixed(2);
}


// And the apir of radio buttons that choose the visitation order.

function chooseVisitationFn(e) {
  visitationFn = this.value;
  if (state === 'running') {
    state = 'switchingVisitationFn';
    clearInterval(timer);
    doRunButton();
  }
}
  

// And finally the radio buttons that choose the acceptance function.

function chooseAcceptanceFn(e) {
  acceptanceFn = this.value;
  if (state === 'running') {
    state = 'switchingAcceptanceFn';
    clearInterval(timer);
    doRunButton();
  }
}
  
// At this point during program load, set the wheels in motion.

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
// temperature. The basis of Figures 12 and 13.

function stepResponse(visitFn, acceptFn, T1, T2, T1_steps, T2_steps, burnin, reps) {
  const prefix = `${visitFn}_${acceptFn}_T${pointless(T1)}_T${pointless(T2)}_steps${T1_steps + T2_steps}_reps${reps} = [`
  itinerary = visitFn;
  acceptanceFn = acceptFn;
  let cachedDrawLattice = drawLattice;
  drawLattice = function(i, j) { ; };
  let data = new Array(T1_steps + T2_steps).fill(0);
  init("scrambled");
  for (let i = 0; i < reps; i++) {
    temperature = T1;
    for (let w = 0; w < burnin; w++) {
      runBatch();
    }
    for (let j = 0; j < T1_steps; j++) {
      runBatch();
      data[j] += Math.abs(magnetization());
    }
    temperature = T2;
    for (let k = 0; k < T2_steps; k++) {
      runBatch();
      data[T1_steps + k] += Math.abs(magnetization());
    }
  }
  drawLattice = cachedDrawLattice;
  let results = data.map(x => (x / reps).toFixed(4));
  outputDiv.textContent = prefix;
  for (let i = 0; i < results.length - 1; i++) {
    outputDiv.textContent += results[i] + ", ";
  }
  outputDiv.textContent += results[results.length - 1] + "]";
  return results;
}

