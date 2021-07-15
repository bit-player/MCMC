
// Explore the influence of boundary-condition rules
// on the Ising model of ferromagnetism.

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

const boundaryButtons = document.querySelectorAll('input[name="boundary"]');
for (let b of boundaryButtons) {
  b.addEventListener('change', chooseBoundary);
}

const algoButtons = document.querySelectorAll('input[name="algorithm"]');
for (let b of algoButtons) {
  b.addEventListener('change', chooseAlgorithm);
}

const highlightCheckbox = document.getElementById("highlight-checkbox");   // "dotted Swiss"
highlightCheckbox.onchange = toggleHighlightState;

const outputDiv = document.getElementById("stats-output");


// Set attributes of the lattice display

const gridSize = 102;              // two extra cells, for the boundary halo
const gridSide = gridSize - 1;
const gridTop = 0;
const gridLeft = 0;
const gridBottom = gridSide;
const gridRight = gridSide;

const activeSize = gridSize - 2;   // these constants define the interior
const activeTop = gridTop + 1;     //   region, exclusive of the boundary
const activeLeft = gridLeft + 1;
const activeBottom = gridBottom - 1; 
const activeRight = gridRight - 1; 

const N = activeSize * activeSize;

const cellSize = theCanvas.width / gridSize;     // note that canvas is 612px (up from 600 in other pgms)
const upColor = "#4b0082";            // indigo
const downColor = "#a297ff";          // mauve
const upBorderColor = "#ad102f";      // crimson
const downBorderColor = "#f16984";    // pink
const zeroBorderColor = "#a0a0a0";    // gray



let highlightFlag = false;                       // boolean for the dotted-Swiss marking

let temperature = Number(theTemperatureSlider.value);


// Global variables for the state machine controlling program execution

let timer;
let boundaryRule = 'wraparound';
let algorithm = 'metropolis';
let state = 'paused';


// bind labels to functions that construct the boundary

const boundaryList = 
  {'wraparound': makeWraparoundBorder,
   'underdog': makeUnderdogBorder,
   'twisted': makeTwistedBorder,
   'up': makeUpBorder,
   'down': makeDownBorder,
   'random': makeRandomBorder,
   'zero': makeZeroBorder,
   'sampled': makeSampledBorder
  }



// build the lattice as an array of arrays of cells

let lattice = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    lattice[i] = new Array(gridSize);
}


// Update the canvas and the two meters. Called after
// every macrostep. Fills in colors only for the active area
// of the lattice, not the borders.

function drawLattice() {
  for (let y = activeTop; y <= activeBottom; y++) {
    for (let x = activeLeft; x <= activeRight; x++) {
      ctx.fillStyle = lattice[x][y] === 1 ? upColor : downColor;
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
  if (highlightFlag) { highlightNeutralSites(); }
  Mreadout.innerHTML = formatReadout(magnetization());
  Rreadout.innerHTML = formatReadout(local_correlations());
}

// Generate the full set of border coordinates; used by all
// the border-drawing functions. We go around the 
// perimeter clockwise, tracing the top, right, bottom,
// and left edges, but skipping over the four corner cells.
// The corners get special treatment as perpetual neutrals,
// because they are not neighbors of any active cell.

function* genBorderXYs() {
  let x = activeLeft, y = gridTop, count = 0;
  const perimeter = (4 * activeSize);
  while (count < perimeter) {
    yield {x, y};
    count++;
    if (count < activeSize) {
      x += 1;
    }
    else if (count < (2 * activeSize)) {
      x = gridRight;
      y += 1;
    }
    else if (count < (3 * activeSize)) {
      y = gridBottom;
      x -= 1;
    }
    else {
      x = gridLeft;
      y -= 1;
    }
  }
}

// Set the four corner cells of the lattice to zero, then
// draw them to the screen.

function makeZeroCorners() {
  ctx.fillStyle = zeroBorderColor;
  for (x of [gridLeft, gridRight]) {
    for (y of [gridTop, gridBottom]) {
      lattice[x][y] = 0;
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
}


// Iterate through all border cells, and paint the screen.

function drawBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    let spinval = lattice[xy.x][xy.y];
    if (spinval === 1) {
      ctx.fillStyle = upBorderColor;
      }
    else if (spinval === -1) {
      ctx.fillStyle = downBorderColor;
    }
    else {
      ctx.fillStyle = zeroBorderColor;
    }
    ctx.fillRect(xy.x * cellSize, xy.y * cellSize, cellSize, cellSize);
  }
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
  for (let x = activeLeft; x <= activeRight; x++) {
    for (let y = activeTop; y <= activeBottom; y++) {
      M += lattice[x][y];
    }
  }
  return M / N;  
}

function local_correlations() {
  let R = 0;
  for (let x = activeLeft; x <= activeRight; x++) {
    for (let y = activeTop; y <= activeBottom; y++) {
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

// The following eight functions handle the initial construcion
// of the border styles. Each one is called when the corresponding
// radio button is activated.


// Neutral option: all border cells have value 0.

function makeZeroBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    lattice[xy.x][xy.y] = 0;
  }
}


// All +1.

function makeUpBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    lattice[xy.x][xy.y] = 1;
  }
}


// All -1.

function makeDownBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    lattice[xy.x][xy.y] = -1;
  }
}


// Random choice between +1 and -1.

function makeRandomBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    let s = coinFlip() ? 1 : -1;
    lattice[xy.x][xy.y] = s;
  }
}


// For each border cell choose a random cell within the active
// area of the lattice. Copy the source cell's state to the
// border cell.

function makeSampledBorder() {
  let XYiter = genBorderXYs();
  for (xy of XYiter) {
    let xSource = 1 + (Math.floor(Math.random() * activeSize));
    let ySource = 1 + (Math.floor(Math.random() * activeSize));
    lattice[xy.x][xy.y] = lattice[xSource][ySource];
  }
}


// The same as makeSampledBorder(), except that we negate the
// value of the selected cell.

function makeUnderdogBorder() {
  let XYiter = genBorderXYs();
  for (xy of XYiter) {
    let xSource = 1 + (Math.floor(Math.random() * activeSize));
    let ySource = 1 + (Math.floor(Math.random() * activeSize));
    lattice[xy.x][xy.y] = -lattice[xSource][ySource];
  }
}


// Each border cell mimics the furthest cell in the same row or column
// on the opposite side of the lattice. 

function makeWraparoundBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    if (xy.y === gridTop) {
      lattice[xy.x][xy.y] = lattice[xy.x][activeBottom];
    }
    if (xy.x === gridRight) {
      lattice[xy.x][xy.y] = lattice[activeLeft][xy.y];
    }
    if (xy.y === gridBottom) {
      lattice[xy.x][xy.y] = lattice[xy.x][activeTop];
    }
    if (xy.x === gridLeft) {
      lattice[xy.x][xy.y] = lattice[activeRight][xy.y];
    }
  }
}


// Each border cell mimics the furthest cell in the same row or column
// on the opposite side of the lattice, but the left-to-right match
// is twisted 180 degrees. 

function makeTwistedBorder() {
  const XYiter = genBorderXYs();
  for (xy of XYiter) {
    if (xy.y === gridTop) {
      lattice[xy.x][xy.y] = lattice[xy.x][activeBottom];
    }
    if (xy.x === gridRight) {
      lattice[xy.x][xy.y] = lattice[activeLeft][gridBottom - xy.y];
    }
    if (xy.y === gridBottom) {
      lattice[xy.x][xy.y] = lattice[xy.x][activeTop];
    }
    if (xy.x === gridLeft) {
      lattice[xy.x][xy.y] = lattice[activeRight][gridBottom - xy.y];
    }
  }
}


// The "make...Border()" functions set up the initial configuration
// for each border style, but updates may be needed as the lattice evolves
// under the action of the Monte Carlo algorithm. In this regard, the
// border styles fall into three classes:
//
//   * the zero, all-up, and all-down borders are static and need 
//     no updating
//
//   * the sampled, underdog, and random borders are updated after
//     each macrostep sweep of the lattice. We can use the 
//     "make...Border()" functions for this purpose.
//
//   * the wraparound and twisted borders must be kept in synch
//     with each incremental change in the lattice, and so they
//     are updated after every microstep.


// Runs after every macrostep. Reconstructs the entire border.

function updateBorderMacrostep() {
  if (boundaryRule === "sampled") {
    makeSampledBorder();
  }
  else if (boundaryRule === "underdog") {
    makeUnderdogBorder();
  }
  else if (boundaryRule === "random") {
    makeRandomBorder();
  }
}


// Runs after every microstep. Updates only the single cell on
// which it is called (and nothing will happen unless that cell
// belongs to the border). Note that we don't have to worry about
// the corner cells because they have value = 0.  

function updateBorderMicrostep(x, y) {
  if (boundaryRule === "wraparound") {
    if (x === activeLeft) {
      lattice[gridRight][y] = lattice[x][y];
    }
    else if (x === activeRight) {
      lattice[gridLeft][y] = lattice[x][y];
    }
    else if (y === activeTop) {
      lattice[x][gridBottom] = lattice[x][y];
    }
    else if (y === activeBottom) {
      lattice[x][gridTop] = lattice[x][y];
    }
  }
  else if (boundaryRule === "twisted") {
    if (x === activeLeft) {
      lattice[gridRight][gridBottom - y] = lattice[x][y];
    }
    else if (x === activeRight) {
      lattice[gridLeft][gridBottom - y] = lattice[x][y];
    }
    else if (y === activeTop) {
      lattice[x][gridBottom] = lattice[x][y];
    }
    else if (y === activeBottom) {
      lattice[x][gridTop] = lattice[x][y];
    }
  }
}


// Draw colored dots on all lattice cells that have exactly two
// friendly and two unfriendly neighbors (neutral, or Swiss, cells).

function highlightNeutralSites() {
  const radius = cellSize / 3;
  const offset = cellSize / 2;
  const tau = 2 * Math.PI;
  for (let x = activeLeft; x <= activeRight; x++) {
    for (let y = activeTop; y <= activeBottom; y++) {
      if (calcDeltaE(x, y) === 0) {
        ctx.fillStyle = (lattice[x][y] === 1) ? "orange" : "green";
        ctx.beginPath();
        ctx.arc(x * cellSize + offset, y * cellSize + offset, radius, 0, tau, true);
        ctx.fill();
      }
    }
  }
}



// Create and display an initial, fully randomized array of cells

function init() {
  let x, y;
  state = 'scrambled';
  for (y = activeTop; y <= activeBottom; y++) {
    for (x = activeLeft; x <= activeRight; x++) {
      lattice[x][y] = coinFlip() ? 1 : -1;
    }
  }
  boundaryList[boundaryRule]();
  makeZeroCorners();
  drawLattice();
  drawBorder();
}



// Given the coordinates of a cell, survey its four neighbors
// and compute the sum of their spins, which will be a value
// in the set {-4, -2, 0, 2, 4}. Multiply by twice the value of the
// central spin, yielding an answer in {-8, -4, 0, 4, 8}. This
// is the net change in system energy if the spin were to flip.

function calcDeltaE(x, y) {                         // assume confined to active grid bounds
  let north, south, east, west, deltaE;
  north = lattice[x][y - 1];
  south = lattice[x][y + 1];
  east  = lattice[x + 1][y];
  west  = lattice[x - 1][y];
  deltaE = 2 * lattice[x][y] * (north + south + east + west);
  return deltaE;
}


// Execute on macrostep sweep of the Matropolis algorithm, using
// the typewriter visitation sequence. Also do the border updates,
// and draw everything to the canvas.

function runMetro() {
  for (let y = activeTop; y <= activeBottom; y++) {
    for (let x = activeLeft; x <= activeRight; x++) {
      let deltaE = calcDeltaE(x, y);
      let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
      if ((deltaE <= 0) || (Math.random() < boltzmann)) {
        lattice[x][y] *= -1;
        updateBorderMicrostep(x, y);
      }
    }
  }
  updateBorderMacrostep();
  drawLattice();
  drawBorder();
}


// Execute on macrostep sweep of the Glauber algorithm, using
// the typewriter visitation sequence. Also do the border updates,
// and draw everything to the canvas.

function runGlauber() {
  for (let i = 0; i < N; i++) {
    let x = 1 + (Math.floor(Math.random() * activeSize));
    let y = 1 + (Math.floor(Math.random() * activeSize));
    let deltaE = calcDeltaE(x, y);
    let boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
    if (Math.random() < (boltzmann / (1 + boltzmann))) {
      lattice[x][y] *= -1;
      updateBorderMicrostep(x, y);
    }
  }
  updateBorderMacrostep();
  drawLattice();
  drawBorder();
}


// Event handler for the paired radio buttons that choose between
// Metropolis and Glauber. If a change is made while the program
// is running, we finish the current macrostep, then restart as
// if the Run button had been pressed.

function chooseAlgorithm(e) {
  algorithm = this.value;
  if (state === 'running') {
    state = 'switchingAlgorithms';
    clearInterval(timer);
    doRunButton();
  }
}


// Event handler for the radio buttons that select a boundary 
// condition. This always calls for redrawing the lattice and the
// border. If a change is made while the program
// is running, we finish the current macrostep, then restart as
// if the Run button had been pressed.

function chooseBoundary(e) {
  boundaryRule = this.value;
  boundaryList[boundaryRule]();
  drawLattice();
  drawBorder();
  if (state === 'running') {
    state = 'switchingBoundaries';
    clearInterval(timer);
    doRunButton();
  }
}
  

// Event handler for the Run/Stop button. 

function doRunButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    this.innerHTML = "Run";
  }
  else {
    state = 'running';
    this.innerHTML = "Stop";
    if (algorithm === "metropolis") {
      timer = setInterval(runMetro, 30);
    }
    else if (algorithm === "glauber") {
      timer = setInterval(runGlauber, 30);
    }
  }
}


// Event handler for single macrostep execution.

function doStepButton(e) {
  if (state === 'running') {
    state = 'paused';
    clearInterval(timer);
    theRunButton.innerHTML = "Run";
  }
  else if (algorithm === "metropolis") {
    runMetro();
  }
  else if (algorithm === "glauber") {
    runGlauber();
  }
}


// Event handler for Reset.

function doResetButton(e) {
  state = 'paused';
  theRunButton.innerHTML = "Run";
  clearInterval(timer);
  outputDiv.textContent = "";
  init();
}


// Event handler for changes to the temperature slider.

function adjustTemperature(e) {
  temperature = Number(this.value);
  temperatureReadout.textContent = temperature.toFixed(2);
}


// Event handler for the checkbox labeled "Dotted Swiss". 

function toggleHighlightState(e) {
  if (highlightFlag) {
    highlightFlag = false;
    drawLattice();
  }
  else {
    highlightFlag = true;
    drawLattice();
  }
}




// On loading, this is where execution begins.

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


// Track the progress of lattice magnetization
// after a sudden change in temperature. This is the code that
// generated Figure 33.

function stepResponse(algo, boundary, T1, T2, T1_steps, T2_steps, burnin, reps) {
  const prefix = `${algo}_${boundary}_T${pointless(T1)}_T${pointless(T2)}_steps${T1_steps + T2_steps}_reps${reps}`;
  let cachedDrawLattice = drawLattice;
  let cachedDrawBorder = drawBorder;
  drawLattice = function() { ; };
  drawBorder = function() { ; };
  let rs = new Array(T1_steps + T2_steps).fill(0);
  let Hs = new Array(T1_steps + T2_steps).fill(0);
  let fn = (algorithm === "metropolis") ? runMetro : runGlauber;
  init();
  for (let i = 0; i < reps; i++) {
    boundaryRule = boundary;
    boundaryList[boundary]()
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
  drawBorder = cachedDrawBorder;
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


// Running tally of how many Swiss cells in indigo and mauve.

function countNeutralSites() {
  let upSum = 0, downSum = 0;
  for (let x = activeLeft; x <= activeRight; x++) {
    for (let y = activeTop; y <= activeBottom; y++) {
      if (calcDeltaE(x, y) === 0) {
        if (lattice[x][y] === 1) {
          upSum++;
        }
        else {
          downSum++;
        }
      }
    }
  }
  console.log(`up neutrals: ${upSum}  down neutrals: ${downSum}`);
}


