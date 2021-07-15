
// A closer look at Monte Carlo algorithms in action, 
// applied to the Ising model of ferromagnetism. Using
// a smaller lattice and a slowed-down pace, we can
// watch an algorithm's path as it progresses from
// site to site, and chooses whether or not to flip
// individual spins. There's also a facility to "doodle"
// patterns on the screen, to see how they fare when
// transformed by the MCMC algorithms.

// MIT License   Copyright (c) 2021 Brian Hayes

// This program accompanies the article "A Month in Monte Carlo"
// at http://bit-player.org/2021/three-months-in-monte-carlo


// Get references to HTML elements; link event handlers

const theViz = document.getElementById("the-viz");
const theCanvas = document.getElementById("the-canvas");
const ctx = theCanvas.getContext("2d");

const Mreadout = document.getElementById("magnetization-readout");
const Rreadout = document.getElementById("correlation-readout");

theCanvas.onclick = doodle;                              // for flipping individual spins
theCanvas.onmouseenter = showCellBoundaries;
theCanvas.onmouseleave = hideCellBoundaries;

const theMicroButton = document.getElementById("micro-button");
theMicroButton.onclick = doMicroButton;
const theMacroButton = document.getElementById("macro-button");
theMacroButton.onclick = doMacroButton;
const theRunButton = document.getElementById("run-button");
theRunButton.onclick = doRunButton;
const theResetButton = document.getElementById("reset-button");
theResetButton.onclick = doResetButton;

const theTemperatureSlider = document.getElementById("temperature-slider");
theTemperatureSlider.onchange = adjustTemperature;
theTemperatureSlider.oninput = adjustTemperature;
const temperatureReadout = document.getElementById("temperature-readout");

const initButtons = document.querySelectorAll('input[name="init-pattern"]');
for (ib of initButtons) {
  ib.addEventListener('change', chooseInitPattern);
}

const visitButtons = document.querySelectorAll('input[name="visit-sequence"]');
for (vb of visitButtons) {
  vb.addEventListener('change', chooseVisitSequence);
}

const acceptButtons = document.querySelectorAll('input[name="acceptance-function"]');
for (ab of acceptButtons) {
  ab.addEventListener('change', chooseAcceptanceFn);
}

const doodleCheckbox = document.getElementById("doodle-checkbox");
doodleCheckbox.onchange = setDoodleMode;

const highlightCheckbox = document.getElementById("highlight-checkbox");
highlightCheckbox.onchange = toggleHighlightState;

const outputDiv = document.getElementById("stats-output");

// Set attributes of the lattice display

const gridSize = 20;
const N = gridSize * gridSize;
const gridEdge = gridSize - 1;
const cellSize = theCanvas.width / gridSize;
const upColor = "#4b0082";
const downColor = "#a297ff";
const upColorFaded = "#423c46";
const downColorFaded = "#c8c7cf";
let temperature = Number(theTemperatureSlider.value);
let stepCount = 0;

let doodleFlag = false;              // true if we're in doodle mode (drawing by hand)
let highlightFlag = false;           // true if showing Dotted Swiss

// One of these will be selected (using the property name)
// to produce a sequence of N pairs of xy coordinates
// for each macrostep of the MCMC algorithm. 

const visitSeqGenerators = {
  'typewriter': typewriterGen,
  'checkerboard': checkerboardGen,
  'diagonal': diagonalGen,
  'permuted': permutedGen,
  'repermuted': re_permutedGen,
  'random': randomGen,
  'simultaneous': simultaneousGen
} 


// And one of these will be chosen to determine whether or
// not each visited spin site should be flipped.

const acceptanceFunctions = {
  'M_rule': MRuleFn,
  'G_rule': GRuleFn,
  'M_star_rule': MstarRuleFn
} 


// Global variables for the state machine controlling program execution

let runTimer;                      // each of these stores a reference
let macrostepTimer;                // to a running timer
let initPattern = 'init-random'    // or all_up, all_down
let visitSeq = 'typewriter'        // or checkerboard, diagonal, permuted, re-permuted, random, simultaneous
let acceptanceFn = 'M_rule'        // or G_rule_, M_star_rule
let state = 'idle';                // just two states: idle and busy

let thePermutation = randomPermutation(N);    // For the permuted and repermuted visitation sequences

let nextXYgen = visitSeqGenerators[visitSeq];      // select one of the sequence generators
let XYiterator = nextXYgen();                      // and get the corresponding iterator

let acceptor = acceptanceFunctions[acceptanceFn];  // similarly, select an acceptance function


// Probability 0.5.

function coinFlip() {
  return Math.random() < 0.5;
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


// Create an array of 1..n, then swap each element in
// turn with a random element. 

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


// Build the lattice as an array of arrays.
// Note there's nothing in these arrays as yet.

let lattice = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    lattice[i] = new Array(gridSize);
}


// A shadow lattice with the same dimensions as the
// main Ising lattice. It will be filled with Booleans
// indicating whether or not a site has already been
// visited.

let visitedArray = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    visitedArray[i] = new Array(gridSize);
}


// Another shadow lattice with the same dimensions as the
// main Ising lattice, for use only in the 'simultaneous'
// visitation sequence. It will be filled with Booleans
// indicating whether or not a site should be flipped at
// the end of the sweep.

let toBeFlippedArray = new Array(gridSize);
for (let i = 0; i < gridSize; i++) {
    toBeFlippedArray[i] = new Array(gridSize);
}


// Color depends on whether the spin is up or down and whether
// to site has been visited before in the current macrostep.

function getSiteColor(x, y) {
  if (lattice[x][y] === 1) {
    if (visitedArray[x][y]) {
      return upColorFaded;        // up, visited
    }
    else {
      return upColor;             // up, not visited
    }
  }
  else {
    if (visitedArray[x][y]) {
      return downColorFaded;      // down, visited
    }
    else {
      return downColor;           // down, not visited
    }
  }
}


// Sweep through all cells of the lattice and paint the canvas
// accordingly. If the highlight flag is set, also add dotted Swiss
// patterns. And the 'simultaneous' visit sequence requires special
// handling: Sites that have not yet flipped but will do so at the
// end of the macrostep need to be marked with a square. Finally
// we update the meters at the top of the display.

function drawLattice() { 
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      ctx.fillStyle = getSiteColor(x, y);
      ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
  if (highlightFlag) { 
    highlightNeutralSites(); 
  }
  if (visitSeq === "simultaneous") {
    markSitesToBeFlipped();
  }
  Mreadout.innerHTML = formatReadout(magnetization());
  Rreadout.innerHTML = formatReadout(local_correlations());
}


// To avoid jitter on the screen, always show a sign prefix,
// and always show three decimal places. (Also use monospace
// font, but that's specified in the CSS.)

function formatReadout(val) {
  sign = (val < 0) ? '&minus;' : '&plus;'
  return sign + Math.abs(val).toFixed(3);
}


// We can use this to reset any of the three lattices -- the
// main one, the visited array, the to-be-flipped array.

function clearLattice(lattice, fillValue) {
  let x, y;
  for (y = 0; y < gridSize; y++) {
    for (x = 0; x < gridSize; x++) {
      lattice[x][y] = fillValue;
    }
  }
}


// Initialize the system upon start or reset.
// Create and display an initial array of cells. The basic
// idea here and in many other procedures is to work with the
// array 'lattice[x][y]' to achieve the configuration we want,
// then call 'drawLattice()' to make it visible.

function init() {
  thePermutation = randomPermutation(N);
  setDoodleMode();
  XYiterator = nextXYgen();
  let x, y;
  if (initPattern === 'init-random') {
    for (y = 0; y < gridSize; y++) {
      for (x = 0; x < gridSize; x++) {
        lattice[x][y] = coinFlip() ? 1 : -1;
      }
    }
  }
  else if (initPattern === 'init-up') {
    for (y = 0; y < gridSize; y++) {
      for (x = 0; x < gridSize; x++) {
        lattice[x][y] = 1;
      }
    }
  }
  else if (initPattern === 'init-down') {
    for (y = 0; y < gridSize; y++) {
      for (x = 0; x < gridSize; x++) {
        lattice[x][y] = -1;
      }
    }
  }
  clearLattice(visitedArray, false);
  clearLattice(toBeFlippedArray, false);
  drawLattice();
  theRunButton.innerHTML = "Run";
}


// Given the coordinates of a cell, survey its four neighbors
// and compute the sum of their spins, which will be a value
// in the set {-4, -2, 0, 2, 4}. Multiply by twice the value of the
// central spin, yielding an answer in {-8, -4, 0, 4, 8}. This
// is the net change in system energy if the spin were to flip.

// This would be a good place to mention that we're using
// periodic boundary conditions: The right edge of the lattice
// wraps around to the left edges, and the bottom wraps around
// to the top.

function calcDeltaE(x, y) {
  let north, south, east, west, deltaE, boltz;
  north = lattice[x][(y > 0) ? y - 1 : gridEdge];
  south = lattice[x][(y < gridEdge) ? y + 1 : 0];
  east  = lattice[(x > 0) ? x - 1 : gridEdge][y];
  west  = lattice[(x < gridEdge) ? x + 1 : 0][y];
  deltaE = 2 * lattice[x][y] * (north + south + east + west);
  return deltaE;
}



// These three functions implement three acceptance rules that
// decide whether a selected spin will be flipped or not. Each rule
// gets the x, y coordinates of a site and returns either true (yes,
// flip the spin), or false (no, leave it as is).

// Note the cutoff in the calculation of the Boltzmann factor: 
// exp(-E/T) overflows IEEE floating point when deltaE = -8 and
// T is small, say 0.1. If the expression returns NaN, we wind
// up with a "stuck pixel" that never flips. Clamping the value
// to the max float avoids this flaw.

// Classic Metropolis. If energetically favorable or neutral, always
// flip. Otherwise, flip with probability proportional to -E/T.

function MRuleFn(x, y) {
  const deltaE = calcDeltaE(x, y);
  const boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
  if (deltaE <= 0) {
    return true;
  }
  else if (Math.random() < boltzmann) {
    return true;
  }
  else {
    return false;
  }
}


// The rule derived from Glauber's 1963 paper. Probability of flipping
// is proportional to a smooth, backwards 'S' function.

function GRuleFn(x, y) {
  const deltaE = calcDeltaE(x, y);
  const boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
  if (Math.random() < (boltzmann / (1 + boltzmann))) {
    return true;
  }
  else {
    return false;
  }
}


// A slight variation on the Metropolis rule. In the neutral case,
// when deltaE === 0, we flip with probability 1/2 rather than 1.

function MstarRuleFn(x, y) {
  const deltaE = calcDeltaE(x, y);
  const boltzmann = Math.min(Math.exp(-deltaE/temperature), Number.MAX_VALUE);
  if (deltaE < 0) {
    return true;
  }
  else if (deltaE === 0) {
    return coinFlip();
  }
  else if (Math.random() < boltzmann) {
    return true;
  }
  else {
    return false;
  }
}



// The visitation sequence has proposed we look at site (x, y).
// We then mark this site as visited, and use one of the three
// acceptor functions described just above to decide whether or
// the spin should be flipped. If the answer is yes, we multiply
// the spin value by -1; otherwise we leave it unchanged. In
// either case we call attention to the current site by drawing
// an outline cursor there.

// Called by doMicroButton, macrostep, and microstep.

function updateSpin(x, y) {  
  visitedArray[x][y] = true;
  let flipSpin = acceptor(x, y);
  if (flipSpin) {
    if (visitSeq === "simultaneous") {
      toBeFlippedArray[x][y] = true;
    }
    else {
      lattice[x][y] *= -1;
    }
  }
  drawLattice();
  drawCursor(x, y);
}
  

// Cursor is a salmon-pink outline drawn along the inner margin
// of the the cell. (Does that make it an inline rather than an outline?)

function drawCursor(x, y) {
  ctx.strokeStyle = "#FF5964";      // hot pink
  ctx.lineWidth = 4;
  ctx.strokeRect((x * cellSize) + 2, (y * cellSize) + 2, cellSize - 4, cellSize - 4);
}


// After various ad hoc solutions in the four earlier programs, I finally
// realized that all the coordinate-generating should be done by function*
// generators. It's a textbook use case.

// The xInit and yInit args are not currently used. They could be used
// to restart the generator from an arbitrary point in the sweep.

function* typewriterGen(xInit=0, yInit=0) {
  for (let y = yInit; y < gridSize; y++) {
    for (let x = xInit; x < gridSize; x++) {
      yield {x, y};
    }
  } 
}

function* simultaneousGen(xInit=0, yInit=0) {     // identical to typewriterGen
  for (let y = yInit; y < gridSize; y++) {
    for (let x = xInit; x < gridSize; x++) {
      yield {x, y};
    }
  }
}

function* checkerboardGen(xInit=0, yInit=0) {
  for (let y = 0; y < gridSize; y++) {
    for (let x = (y % 2); x < gridSize; x += 2) {
      yield {x, y};
    }
  }
  for (let y = 0; y < gridSize; y++) {
    for (let x = ((y + 1) % 2); x < gridSize; x += 2) {
      yield {x, y};
    }
  }
}

function* diagonalGen(xInit=0, yInit=0) {
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

function* permutedGen() {
  for (let r of thePermutation) {
    let x = Math.floor(r / gridSize);
    let y = Math.floor(r % gridSize);
    yield {x, y};
  }
}

function* re_permutedGen() {
  const p = randomPermutation(N);
  for (let r of p) {
    let x = Math.floor(r / gridSize);
    let y = Math.floor(r % gridSize);
    yield {x, y};
  }
}

function* randomGen() {
  for (let i = 0; i < N; i++) {
    let x = Math.floor(Math.random() * gridSize);
    let y = Math.floor(Math.random() * gridSize);
    yield {x, y};
  }
}


// During the lattice sweep under the 'simultaneous' visit sequence, we
// don't actually change any spins, but we need to show which ones
// will be changed when the sweep concludes. This draws a small,
// contrasting square in the middle of each such cell. 

function markSitesToBeFlipped() {
  const markSize = 10;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      const rectLeft = (x * cellSize) + (cellSize - markSize) / 2;
      const rectTop = (y * cellSize) + (cellSize - markSize) / 2;
      if (toBeFlippedArray[x][y]) {
        ctx.fillStyle = (lattice[x][y] === 1) ? downColorFaded : upColorFaded;
        ctx.fillRect(rectLeft, rectTop, markSize, markSize);
      }
    }
  }
}


// Sweep through the lattice once. Called by doMacroButton.

function macrostep() {
  let nextCoords = XYiterator.next();
  if (nextCoords.done) {
    state = 'idle';
    clearInterval(macrostepTimer);
    XYiterator = nextXYgen();
    clearLattice(visitedArray, false);
    if (visitSeq === "simultaneous") {
      simultaneousUpdate();
    }
  }
  else {
    let {x: x, y: y} = nextCoords.value;
    updateSpin(x, y);
  }
}


// The 'simultaneous' routine marks cells to be flipped but does not
// make any changes to the lattice. This is called at the end of each
// sweep to flip all the spins flagged in the 'toBeFlippedArray'.

function simultaneousUpdate() {
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      if (toBeFlippedArray[x][y]) {
        lattice[x][y] *= -1;
      }
    }
  }
  clearLattice(toBeFlippedArray, false);
  drawLattice();
}


// Run continuous sweeps of the chosen visitation sequence. Called at
// 10-millisecond intervals by doRunButton.

function runForever() {
  let nextCoords = XYiterator.next();
  if (nextCoords.done) {
    XYiterator = nextXYgen();
    nextCoords = XYiterator.next();
    clearLattice(visitedArray, false);
    if (visitSeq === "simultaneous") {
      simultaneousUpdate();
    }
  }
  let {x: x, y: y} = nextCoords.value;
  updateSpin(x, y);
}


// Event handler for the Microstep button. Moves the cursor to
// the next cell in the XYiterator sequence, restarting that
// sequence in the case that done===true. Then calls updateSpin(x, y)
// to do the rest.

function doMicroButton(e) {
  if (state === 'busy') {
    clearInterval(runTimer);
    clearInterval(macrostepTimer);
    theRunButton.innerHTML = "Run";
    state = 'idle';
  }
  else {
    let nextCoords = XYiterator.next();
    if (nextCoords.done) {
      clearLattice(visitedArray, false);
      clearLattice(toBeFlippedArray, false);
      XYiterator = nextXYgen();
      nextCoords = XYiterator.next();
    }
    let {x: x, y: y} = nextCoords.value;
    drawLattice();                        // to erase previous highlight
    updateSpin(x, y);
  }
}


// Event handler for the MacroStep button. If either macrostep or
// runForever is currently running, shuts them down and waits for
// further instructions. On the other hand, if we're at idle,
// call 'macrostep' at 50-millisecond intervals. 

function doMacroButton(e) {
  if (state === 'busy') {
    state = 'idle';
    clearInterval(runTimer);
    clearInterval(macrostepTimer);
    theRunButton.innerHTML = "Run";
  }
  else {
    state = 'busy';
    macrostepTimer = setInterval(macrostep, 50);
  }
}


// Event handler for the Run/Stop button. If either macrostep or
// runForever is currently running, shuts them down and waits for
// further instructions. On the other hand, if we're at idle,
// call 'runForever' at 10-millisecond intervals. 

function doRunButton(e) {
  if (state === 'busy') {
    state = 'idle';
    clearInterval(runTimer);
    clearInterval(macrostepTimer);
    this.innerHTML = "Run";
  }
  else {
    state = 'busy';
    this.innerHTML = "Stop";
    runTimer = setInterval(runForever, 10);
  }
}


// Event handler for the Reset button. Shuts down timers and
// calls init().

function doResetButton(e) {
  state = 'idle';
  clearInterval(runTimer);
  clearInterval(macrostepTimer);
  init();
}


// Event handler for the temperature slider.

function adjustTemperature(e) {
  temperature = Number(this.value);
  temperatureReadout.textContent = temperature.toFixed(2);
}


// Handler for change event in the radio button group
// for visitation sequence. Binds a new generator to
// the variable nextXYgen, then sets this function as the
// new value of XYiterator. Thus on the next call to
// Micro, Macro, or Run, the next lattice site will come
// from this new sequence, even if we've interrupted another
// sequence in mid-sweep.

function chooseVisitSequence(e) {
  visitSeq = this.value;
  nextXYgen = visitSeqGenerators[visitSeq];
  XYiterator = nextXYgen();
}


// Event handler for group of radio buttons.
// Changing acceptance function takes effect immediately.

function chooseAcceptanceFn(e) {
  acceptanceFn = this.value;
  acceptor = acceptanceFunctions[acceptanceFn];
}


// Event handler for group of radio buttons.
// Choosing a new init pattern shuts down everything and
// calls init() with the new pattern in effect.

function chooseInitPattern(e) {
  initPattern = this.value;
  if (state === 'idle') {
    clearInterval(macrostepTimer);
    clearInterval(runTimer);
    init();
  }
}


// Event handler for the checkbox that turns dotted-Swissing
// on and off.

function toggleHighlightState(e) {
  highlightFlag = !highlightFlag
  drawLattice();
}


// The function called by drawLattice to do the Swiss dotting,
// but only if the highlightFlag is lit. Puts orange dots in
// indigo cells and green dots in mauve cells.

function highlightNeutralSites() {
  const radius = cellSize / 3;
  const offset = cellSize / 2;
  const tau = 2 * Math.PI;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      if (calcDeltaE(x, y) === 0) {
        ctx.fillStyle = (lattice[x][y] === 1) ? "orange" : "green";
        ctx.beginPath();
        ctx.arc(x * cellSize + offset, y * cellSize + offset, radius, 0, tau, true);
        ctx.fill();
      }
    }
  }
}

// Event handler for the checkbox that activates doodle mode.
// Doodle mode allows you to set the state of individual spins by
// clicking on the canvas. 

function setDoodleMode(e) {
  if (doodleCheckbox.checked) {
    doodleFlag = true;
  }
  else {
    doodleFlag = false;
  }
}


// Event handler for clicks within the canvas. If we're in
// doodle mode, we figure out where the click occurred and
// toggle the state of the underlying cell. (If we're not
// in doodle mode, do nothing.)

function doodle(e) {
  if (doodleFlag) {
    const rect = theCanvas.getBoundingClientRect()
    const pixelX = event.clientX - rect.left
    const pixelY = event.clientY - rect.top
    const cellX = Math.floor(pixelX / cellSize);
    const cellY = Math.floor(pixelY / cellSize);
    lattice[cellX][cellY] *= -1;
    drawLattice();
    showCellBoundaries();
  }
}


// Event handler for mouseenter events on the canvas. Shows
// a grid of light gray lines outlining the cells. An aid for
// doodle mode, although I've decided to draw the grid even
// when doodling is not active.

function showCellBoundaries(e) {
  ctx.strokeStyle = "#ccc";
  ctx.lineWidth = 1;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      ctx.strokeRect(x * cellSize, y * cellSize, cellSize, cellSize);
    }
  }
}


// Event handler for mouseleave events on the canvas. Erase
// the grid.

function hideCellBoundaries(e) {
  drawLattice();
}
  

// During program loading, at this point we call the init()
// function. Once it returns we're set to go.

init();

// ******* END OF INTERACTIVE SIMULATION CODE *******
// ******* All functions below this line must *******
// ******* be invoked from the console.       *******


// Sometimes we want to save an image of some state of the lattice.
// Executing this function will add a link that can be clicked to
// save the current state as a PNG file.

function showSaveLink() {
  vizDiv = document.getElementById("the-viz");
  const saveLink = document.createElement("a");
  saveLink.href = "dummy";
  saveLink.download = "ising_lattice.png";
  saveLink.innerHTML = "Save Canvas";
  saveLink.onclick = saveCanvas;
  vizDiv.appendChild(saveLink);
}


// The event handler invoked by the save link.

function saveCanvas(e) {
  canvasData = theCanvas.toDataURL();
  this.href = canvasData; 
}


// Running tally of how many Swiss cells in indigo and mauve.

function countNeutralSites() {
  let upSum = 0, downSum = 0;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
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


// Doodle mode lets you draw pixel by pixel, but the clicking gets
// tiresome. Herewith canned routines for a few simple shapes.
// These are designed to draw a shape made up of indigo cells on
// a mauve background.

function makeBlock(side) {
  const offset = Math.floor((gridSize - side) / 2);
  for (let y = offset; y < offset + side; y++) {
    for (let x = offset; x < offset + side; x++) {
      lattice[x][y] = 1;
    }
  }
  drawLattice();
}

function makeRiver(width) {
  const offset = Math.floor((gridSize - width) / 2);
  for (let y = offset; y < (offset + width); y++) {
    for (let x = 0; x < gridSize; x++) {
      lattice[x][y] = 1;
    }
  }
  drawLattice();
}

function makeSinusoidalRiver(wavenumber, riverWidth, meanderAmplitude) {

  const tau = 2 * Math.PI; 
  const sine = function (x) {
    return Math.floor(meanderAmplitude * Math.sin((x / gridSize) * tau * wavenumber));
  }
  
  clearLattice(lattice, -1);
  
  for (let x = 0; x < gridSize; x++) {
    let sinx = sine(x);
    let miny = (gridSize / 2) - (riverWidth / 2) + sinx;
    let maxy = (gridSize / 2) + (riverWidth / 2) + sinx;
    for (let y = miny; y <= maxy; y++) {
      lattice[x][y] = 1;
    }
  }
  drawLattice();
}

function makeCircle(radius) {

  const center = gridSize / 2;
  
  function distanceToCenter(x, y) {
    const dx = x - center;
    const dy = y - center;
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


// Tally cells in four categories.

function tallyUpsDownsNeutrals() {
  let upCells = 0;
  let downCells = 0;
  let upNeutrals = 0;
  let downNeutrals = 0;
  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      if (lattice[x][y] === 1) {
        upCells++
        if (calcDeltaE(x, y) === 0) {
          upNeutrals++
        }
      }
      else {
        downCells++
        if (calcDeltaE(x, y) === 0) {
          downNeutrals++
        }
      }
    }
  }
  return {uc: upCells, dc: downCells, un: upNeutrals, dn: downNeutrals};
}


// Code for producing Figures 26, 27, and 30.

function traceNeutralCounts(stepLimit, filenamePrefix) {
  
  function tally() {
    let upCells = 0;
    let downCells = 0;
    let upNeutrals = 0;
    let downNeutrals = 0;
    for (let y = 0; y < gridSize; y++) {
      for (let x = 0; x < gridSize; x++) {
        if (lattice[x][y] === 1) {
          upCells++
          if (calcDeltaE(x, y) === 0) {
            upNeutrals++
          }
        }
        else {
          downCells++
          if (calcDeltaE(x, y) === 0) {
            downNeutrals++
          }
        }
      }
    }
    indexList.push(stepCount);
    upCellsList.push(upCells);
    downCellsList.push(downCells);
    upNeutralsList.push(upNeutrals);
    downNeutralsList.push(downNeutrals);
    return Math.min(upCells, downCells);
  }
  
  let stepCount = 0;
  const indexList = []
  const upCellsList = [];
  const downCellsList = [];
  const upNeutralsList = [];
  const downNeutralsList = [];
  lesserCellCount = tally();
  while ((stepCount < stepLimit) && (lesserCellCount > 0)) {
    stepCount++;
    let nextCoords = XYiterator.next();
    if (nextCoords.done) {
      XYiterator = nextXYgen();
      nextCoords = XYiterator.next();
    }
    let {x: x, y: y} = nextCoords.value;
    let flipSpin = acceptor(x, y);
    if (flipSpin) {
      lattice[x][y] *= -1;
      lesserCellCount = tally();
    }
  }
  const neutralDiffList = [];
  for (let i = 0; i < upCellsList.length; i++) {
    neutralDiffList[i] = upNeutralsList[i] - downNeutralsList[i];
  }
  
  outputDiv.textContent = filenamePrefix + "_indices = [";
  outputDiv.textContent += indexList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_upCells = [";
  outputDiv.textContent += upCellsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_downCells = [";
  outputDiv.textContent += downCellsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_upNeutrals = [";
  outputDiv.textContent += upNeutralsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_downNeutrals = [";
  outputDiv.textContent += downNeutralsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_neutralDiffs = [";
  outputDiv.textContent += neutralDiffList.join(", ");
  outputDiv.textContent += "] ";
}




// Late in the game, I finally realized I could save a great deal of
// computing time by running the MCMC lattice sweeps in a simple
// while loop, rather than calling the functions tied to the graphic
// interface.

function statsNeutralCounts(initFn, stepLimit, filenamePrefix, reps) {
  
  function tally() {
    let upCells = 0;
    let downCells = 0;
    let upNeutrals = 0;
    let downNeutrals = 0;
    for (let y = 0; y < gridSize; y++) {
      for (let x = 0; x < gridSize; x++) {
        if (lattice[x][y] === 1) {
          upCells++
          if (calcDeltaE(x, y) === 0) {
            upNeutrals++
          }
        }
        else {
          downCells++
          if (calcDeltaE(x, y) === 0) {
            downNeutrals++
          }
        }
      }
    }
    upCellsList[stepCount] += upCells;
    downCellsList[stepCount] += downCells;
    upNeutralsList[stepCount] += upNeutrals;
    downNeutralsList[stepCount] += downNeutrals;
  }
  
  function macroQuickstep() {
    let XYiterator = nextXYgen();
    let nextCoords = XYiterator.next();
    while (!nextCoords.done) {
      let {x: x, y: y} = nextCoords.value;
      let flipSpin = acceptor(x, y);
      if (flipSpin) {
        lattice[x][y] *= -1;
      }
      nextCoords = XYiterator.next();
    }
  }
  
  let upCellsList = new Array(stepLimit).fill(0);
  let downCellsList = new Array(stepLimit).fill(0);
  let upNeutralsList = new Array(stepLimit).fill(0);
  let downNeutralsList = new Array(stepLimit).fill(0);
  
  for (i=0; i<reps; i++) {
    // clearLattice(lattice, -1);
    initFn();
    for (stepCount=0; stepCount<stepLimit; stepCount++) {
      tally();
      macroQuickstep();
    }
  }  
  // console.log(upCellsList, downCellsList, upNeutralsList, downNeutralsList);
  let neutralDiffList = [];
  for (let i = 0; i < upCellsList.length; i++) {
    neutralDiffList[i] = upNeutralsList[i] - downNeutralsList[i];
  }
  
  upCellsList = upCellsList.map(x => (x / reps).toFixed(3));
  downCellsList = downCellsList.map(x => (x / reps).toFixed(3));
  upNeutralsList = upNeutralsList.map(x => (x / reps).toFixed(3));
  downNeutralsList = downNeutralsList.map(x => (x / reps).toFixed(3));
  neutralDiffList = neutralDiffList.map(x => (x / reps).toFixed(3));
  
  outputDiv.textContent = "";
  
  outputDiv.textContent += filenamePrefix + "_upCells = [";
  outputDiv.textContent += upCellsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_downCells = [";
  outputDiv.textContent += downCellsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_upNeutrals = [";
  outputDiv.textContent += upNeutralsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_downNeutrals = [";
  outputDiv.textContent += downNeutralsList.join(", ");
  outputDiv.textContent += "] ";
  
  outputDiv.textContent += filenamePrefix + "_neutralDiffs = [";
  outputDiv.textContent += neutralDiffList.join(", ");
  outputDiv.textContent += "] ";
}


// The code that generated Figure 32. (But with lattice size
// set to 100 x 100.)

function lifetimes(initFn, filename, stepLimit, reps) {
  
  function macroQuickstep() {
    let XYiterator = nextXYgen();
    let nextCoords = XYiterator.next();
    while (!nextCoords.done) {
      let {x: x, y: y} = nextCoords.value;
      let flipSpin = acceptor(x, y);
      if (flipSpin) {
        lattice[x][y] *= -1;
      }
      nextCoords = XYiterator.next();
    }
  }
  
  function done() {
    let upCount = 0;
    for (let y = 0; y < gridSize; y++) {
      for (let x = 0; x < gridSize; x++) {
        if (lattice[x][y] === 1) {
          upCount++;
        }
      }
    }
    return ((upCount < gridSize) || (upCount > (N - gridSize)));
  }
  
  outputDiv.textContent = "";
  const datalist = [];
  for (i = 0; i < reps; i++) {
    console.log(i);
    clearLattice(lattice, -1)
    initFn();
    stepCount = 0;
    while (!done() && stepCount < stepLimit) {
      macroQuickstep();
      stepCount++
    }
    datalist.push(stepCount);
  }
  outputDiv.textContent += `${filename} = [`;
  outputDiv.textContent += datalist.join(", ");
  outputDiv.textContent += "] ";
}


