
// Simulate Glauber dynamics and the Metropolis algorithm
// in the Ising model of ferromagnetism.

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

const radioButtons = document.querySelectorAll('input[name="algorithm"]');
for (rb of radioButtons) {
  rb.addEventListener('change', chooseAlgorithm);
}

const outputDiv = document.getElementById("stats-output");

// Set attributes of the lattice display

const gridSize = 100;            // cells per row and per column
const N = gridSize * gridSize;   // = 10,000
const gridEdge = gridSize - 1;   // last index in 0..99
const cellSize = theCanvas.width / gridSize;
const upColor = "#4b0082";           // indigo
const downColor = "#a297ff";         // call it mauve, though it's really light blue
let temperature = Number(theTemperatureSlider.value);


// Global variables for the state machine controlling program execution

let timer;
let algorithm = 'metro'   // other option is 'glauber'
let state = 'paused';


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

const square = x => x * x;


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
  let north, south, east, west, deltaE, boltz;
  north = lattice[x][(y > 0) ? y - 1 : gridEdge];
  south = lattice[x][(y < gridEdge) ? y + 1 : 0];
  east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
  west  = lattice[(x < gridEdge) ? x + 1 : 0][y];
  return 2 * lattice[x][y] * (north + south + east + west);
}


// Make a single sweep through the entire lattice, updating
// spins according to the Metropolis algorithm, then redraw
// the lattice. 

// NOTE ON BOLTZMANN WEIGHT: After calculating the weight as
// exp(dE/T), why the funny business about Math.min and 
// MAX_VALUE? Because if dE has a large negative value,
// and T is very small, exp(dE/T) can cause floating-point
// overflow. For example, exp(8/0.01) is > 10^347. JavaScript
// barfs on anything greater than 10^308.

function runMetro() {
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      let deltaE = calcDeltaE(x, y);
      let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
      if ((deltaE <= 0) || (Math.random() < boltzmann)) {     // this is the M-rule
        lattice[x][y] *= -1;
      }
    }
  }
  drawLattice();
}


// The corresponding procedure for Glauber dynamics. Generates
// N pairs of random xy coordinates.

function runGlauber() {
  for (let i = 0; i < N; i++) {
    let x = Math.floor(Math.random() * gridSize);
    let y = Math.floor(Math.random() * gridSize);
    let deltaE = calcDeltaE(x, y);
    let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
    if (Math.random() < (boltzmann / (1 + boltzmann))) {      // the G-rule
      lattice[x][y] *= -1;
    }
  }
  drawLattice();
}


// Event handler for presses of the Run/Stop button. This
// is a toggle switch: It starts the simulation if the
// program is idle, but if the simulation is already running,
// the button stops it.

// The program runs in a timer loop, using setInterval to
// start a new instance every 30 milliseconds.

function doRunButton(e) {
  if (state !== 'running') {
    state = 'running';
    this.innerHTML = "Stop";
    if (algorithm === 'metro') {
      timer = setInterval(runMetro, 30);
    }
    else if (algorithm === 'glauber') {
      timer = setInterval(runGlauber, 30);
    }
  }
  else {
    state = 'paused';
    clearInterval(timer);
    this.innerHTML = "Run";
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
    if (algorithm === 'metro') {
      runMetro();
    }
    else if (algorithm === 'glauber') {
      runGlauber();
    }
  }
}


// Event handler for the Reset button. Clears the timer
// to stop the animation, then re-inits the program.
// Also clears the outputDiv, about which more below,
// under console-only functions.

function doResetButton(e) {
  state = 'paused';
  theRunButton.innerHTML = "Run";
  clearInterval(timer);
  outputDiv.textContent = "";
  init();
}


// Event handler for the temperature slider. Update the
// 'temperature' global, and the readout for the slider.

function adjustTemperature(e) {
  temperature = Number(this.value);
  temperatureReadout.textContent = temperature.toFixed(2);
}


// Event handler for the group of two radio buttons that
// select either Matropolis or Glauber algorithms. If the
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
// function. Once it returns we're set yo go.

init();


// ******* END OF INTERACTIVE SIMULATION CODE *******
// ******* All functions below this line must *******
// ******* be invoked from the console.       *******




// Sum up the energy in the 2N bonds between adjacent
// spins in the lattice. Each pair with the same orientation
// counts as -1, and each opposite pair is +1.

function calc_lattice_energy() {
  let s = 0;
  let east, south;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      here = lattice[x][y];
      east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
      south = lattice[x][(y < gridEdge) ? y + 1 : 0];
      if (here === east) {
        s -= 1;
      }
      else {
        s += 1
      }
      if (here === south) {
        s -= 1;
      }
      else {
        s += 1
      }
    }
  }
  return s;
}


// Hmm. What does this do, and what did I use it for?

function random_state_energy_distribution(reps) {
  let E, absE;
  const pluses = new Array(200).fill(0);
  const minuses = new Array(200).fill(0);
  let zeroCount = 0;
  for (i = 0; i < reps; i++) {
    init();
    E = calc_lattice_energy();
    absE = Math.floor(Math.abs(E));
    if (E > 0) {
      if (pluses[absE]) {
        pluses[absE] += 1;
      }
      else {
        pluses[absE] = 1;
      }
    }
    if (E < 0) {
      if (minuses[absE]) {
        minuses[absE] += 1;
      }
      else {
        minuses[absE] = 1;
      }
    }
    if (E === 0) {
      zeroCount++
    }
  }
  return {minuses, zeroCount, pluses};
}


// The 'dotted Swiss' function.

function highlightNeutralSites() {
  const radius = cellSize / 3;
  const offset = cellSize / 2;
  const twopi = 2 * Math.PI;
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      if (calcDeltaE(x, y) === 0) {
        ctx.fillStyle = (lattice[x][y] === 1) ? "orange" : "green";
        ctx.beginPath();
        ctx.arc(x * cellSize + offset, y * cellSize + offset, radius, 0, twopi, true);
        ctx.fill();
      }
    }
  }
}


// See Figure 11 in the article. Recolor the lattice sites based
// on the number of friendly neighbors. 

function recolorByLocalCorrelation() {
  const palette = ["#810f7c", "#8856a7", "#8c96c6", "#b3cde3", "#edf8fb"];
  for (let x = 0; x < gridSize; x++) {
    for (let y = 0; y < gridSize; y++) {
      let paletteIndex = (calcDeltaE(x, y) / 4) + 2;
      let color = palette[paletteIndex];
      ctx.fillStyle = color;
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
}




// Convert a number to a string, replacing
// any decimal point with "p", to make it suitable
// as a component of a file name.

function pointless(num) {
  return String(num).replace(".", "p");
}


// The function used to produce data for the table titled
// "Magnitude of Magnetization."

// Output goes to a div included in the DOM after the canvas
// and its control panel. That div is originally empty and
// hence invisible. This is a cheap trick, adopted mainly
// because output to the console is hard to copy and paste.

function equilibriumProperties(algo, initState, T, burnin, steps, reps) {
  const prefix = `${algo}_${initState}_${pointless(T)}_${burnin}_${steps}_${reps}`
  let fn = (algo === "metro") ? runMetro : runGlauber;
  let Rs = new Array(steps).fill(0);   // local correlation values
  let Ms = new Array(steps).fill(0);   // magnetization values
  let Es = new Array(steps).fill(0);   // total lattice energy values
  
  let cachedDrawLattice = drawLattice;   // nullify drawLattice fn, so we
  drawLattice = function(i, j) { ; };    // don't waste time updating canvas
  
  temperature = T;
  for (let i = 0; i < reps; i++) {       // The main loop. For each rep,
    init(initState);                     // initialise the lattice,
    for (let w = 0; w < burnin; w++) {   // throw away 'burnin' macrosteps
      fn();
    }
    for (let j = 0; j < steps; j++) {    // collect data for 'steps' steps
      fn();
      Es[j] += calc_lattice_energy();
      Rs[j] += local_correlations();
      Ms[j] += Math.abs(magnetization());
    }
  }
  drawLattice = cachedDrawLattice;    // restore original canvas function
  
  let E_results = Es.map(x => (x / reps));   // reduce to [0,1] range
  let R_results = Rs.map(x => (x / reps));
  let M_results = Ms.map(x => (x / reps));  
  
  // Write the list of lattice energies to the output div, after
  // first sorting them into descending order. I don't actually
  // remember why I wanted to do this -- maybe just to
  // check that the range of values is not too wide, that our
  // average values are also typical values.
  
  let E_list = "<p>" + prefix + "_E = [";
  E_results.sort((a,b) => a - b);
  for (let i = 0; i < E_results.length - 1; i++) {
    E_list += E_results[i].toString() + ", ";
  }
  E_list += E_results[E_results.length - 1].toString() + "]</p>";  // no comma at end of list
  outputDiv.innerHTML = E_list;
 
  // now do the stats
  const count = reps * steps;
  const R_mean = R_results.reduce((acc, x) => acc + x, 0) / steps;
  const R_deviates = R_results.map(x => square(x - R_mean));
  const R_sdev = Math.sqrt(R_deviates.reduce((sum, x) => sum + x) / (steps - 1));
  const M_mean = M_results.reduce((acc, x) => acc + x, 0) / steps;
  const M_deviates = M_results.map(x => square(x - M_mean));
  const M_sdev = Math.sqrt(M_deviates.reduce((sum, x) => sum + x) / (steps - 1));
  
  outputDiv.innerHTML += `<p>R = ${R_mean} ± ${R_sdev}</p>`;
  outputDiv.innerHTML += `<p>M = ${M_mean} ± ${M_sdev}</p>`;
  
  return {R_mean, R_sdev, M_mean, M_sdev};  // goes to console
}


// The first version of a routine that appears in later programs
// as well, meant to track the progress of lattice magnetization
// after a sudden change in temperature. This the code that
// generated Figure 8.

function temperatureStep(algo, initState, T1, T2, T1_steps, T2_steps, glom, burnin, reps) {
  const prefix = `${algo}_T${pointless(T1)}_T${pointless(T2)}_steps${T1_steps + T2_steps}_reps${reps}`;  // for file name
  let cachedDrawLattice = drawLattice;
  drawLattice = function(i, j) { ; };                 // no drawing to canvas during this routine
  let Rs = new Array(T1_steps + T2_steps).fill(0);
  let Ms = new Array(T1_steps + T2_steps).fill(0);
  let fn = (algo === "metro") ? runMetro : runGlauber;
  for (let i = 0; i < reps; i++) {                       // run routine 'reps' times
    init(initState);                                     // initialize lattice
    temperature = T1;                                    // starting temperature
    for (let w = 0; w < burnin; w++) {
      fn();                                              // throw away burnin data
    }
    for (let j = 0; j < T1_steps; j++) {                 // gather data at T1
      fn();
      Rs[j] += local_correlations();
      Ms[j] += Math.abs(magnetization());
    }
    temperature = T2;                                    // set temperature to T2
    for (let k = 0; k < T2_steps; k++) {                 // gather data at T2
      fn();
      Rs[T1_steps + k] += local_correlations();
      Ms[T1_steps + k] += Math.abs(magnetization());
    }
  }
  drawLattice = cachedDrawLattice;                       // restore cached fn def
  
  let R_results = Rs.map(x => (x / reps));               // reduce to [0,1] range
  let M_results = Ms.map(x => (x / reps));
  let R_short = decimate_array(R_results, glom);         // very long runs produce too much data; consolidate
  let M_short = decimate_array(M_results, glom);
  
  let R_str = "<p>" + prefix + "_R = [";                 // write outputs to outputDiv
  for (let i = 0; i < R_short.length - 1; i++) {
    R_str += R_short[i].toFixed(4) + ", ";
  }
  R_str += R_short[R_short.length - 1].toFixed(4) + "]</p>";  // no comma at end of list
  outputDiv.innerHTML = R-str;
 
  let M_str = "<p>" + prefix + "_M = [";
  for (let i = 0; i < M_short.length - 1; i++) {
    M_str += M_short[i].toFixed(4) + ", ";
  }
  M_str += M_short[M_short.length - 1].toFixed(4) + "]</p>";  // no comma at end of list
  outputDiv.innerHTML = M-str;
}


// Let g = glom_size. Take every consecutive block of g
// elements and calculate its average value. Construct a
// new array as the sequence of averages.

// (What if arr.length is not a multiple of glom_size?
// You'll get some NaN values at the end. Sorry. Lazy.)

function decimate_array(arr, glom_size) {
  const shortArray = [];
  for (outer = 0; outer < arr.length; outer += glom_size) {
    let sum = 0;
    for (inner = 0; inner < glom_size; inner++) {
      sum += arr[outer + inner];
    }
    shortArray.push(sum / glom_size);
  }
  return shortArray;
}


// Sometimes we want to save an image of some state of the lattice.
// Executing this function will add a link that can be clicked to
// save the current state as a PNG file.

function showSaveLink() {
  const saveLink = document.createElement("a");
  saveLink.href = "dummy";
  saveLink.download = "ising_lattice.png";
  saveLink.innerHTML = "Save Canvas";
  saveLink.onclick = saveCanvas;
  theViz.appendChild(saveLink);
}

// The event handler invoked by the save link.

function saveCanvas(e) {
  canvasData = theCanvas.toDataURL();
  this.href = canvasData; 
}



// Code for generating random arrays of spins, then seeing
// whether any among them might be expected to turn up at the
// given temperature. 

function hitOrMiss(T, reps) {
  let cachedDrawLattice = drawLattice;   // nullify drawLattice fn, so we
  drawLattice = function(i, j) { ; };    // don't waste time updating canvas
  const prefix = `hit_or_miss_${pointless(T)}_${reps}`
  let weights = [];
  let ETs = [];
  let probs = [];
  let maxWeight = 0;
  let minWeight = Infinity;
  let sumOfWeights = 0;
  for (i=0; i<reps; i++) {
    init('scrambled');
    let E = calc_lattice_energy();
    let EoverT = E/T;
    ETs[i] = EoverT;
    wt = Math.exp(-EoverT);
    weights[i] = wt;
    maxWeight = Math.max(wt, maxWeight);
    minWeight = Math.min(wt, minWeight);
    sumOfWeights += wt
  }
  let meanWeight = sumOfWeights / reps;
  for (j=0; j<reps; j++) {
    probs[j] = weights[j] / sumOfWeights;
  }
  drawLattice = cachedDrawLattice;                       // restore cached fn def

  ETs.sort((a,b) => a - b);
  outputDiv.textContent += prefix + "_E = [";
  for (let i = 0; i < ETs.length - 1; i++) {
    outputDiv.textContent += ETs[i].toString() + ", ";
  }
  outputDiv.textContent += ETs[ETs.length - 1].toString() + "]";
 
  return {maxWeight, meanWeight, minWeight, ETs, weights, probs};
}




// Generate data for the long-tailed histogram of convergence
// times in Figure 18.

function convergence_time_histogram(algo, threshold, reps) {
  outputDiv.textContent = "";
  const prefix = `convergencetime_${algo}_${pointless(threshold)}_${reps}`
  const histo_data = [];
  let fn = (algo === "metro") ? runMetro : runGlauber;
  let cachedDrawLattice = drawLattice;
  drawLattice = function(i, j) { ; };

  for (i=0; i<reps; i++) {
    init('scrambled');
    temperature = 10;
    for (let w = 0; w < 200; w++) {
      fn();
    }
    temperature = 2;
    let M = magnetization();
    stepcount = 0;
    while ((Math.abs(M) < threshold) && (stepcount < 10000)) {
      // console.log(stepcount, M);
      fn();
      stepcount++
      M = magnetization();
    }
    histo_data.push(stepcount);    
  }
  drawLattice = cachedDrawLattice;
  outputDiv.textContent = prefix + " = [";
  for (let i = 0; i < histo_data.length - 1; i++) {
    outputDiv.textContent += histo_data[i].toString() + ", ";
  }
  outputDiv.textContent += histo_data[histo_data.length - 1].toString() + "]";
}


// This routine was useful for creating the circular blobs in
// Figures 19 and 20.

function makeCircle(radius, center_x, center_y) {

  function distanceToCenter(x, y) {
    const dx = Math.min(x - center_x, center_x - (x - gridSize));
    const dy = Math.min(y - center_y, center_y - (y - gridSize));
    // console.log(dx, dy);
    return Math.sqrt((dx * dx) + (dy * dy));
  }
  
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      if (distanceToCenter(x, y) <= radius) {
        lattice[x][y] = 1;
      }
    }
  }
  drawLattice();
}

