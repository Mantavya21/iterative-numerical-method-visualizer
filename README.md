**Overview:**
The Iterative Numerical Methods Visualizer is a client-side web application designed to demonstrate the geometric behavior and convergence patterns of common root-finding and interpolation algorithms.

Instead of acting as a simple calculator that outputs a final answer, this tool allows users to step through computational mathematics frame-by-frame. It is built to help visualize edge cases, algorithm breakdown points, and the mathematical logic behind iterative numerical analysis.

**Features:**
1) Dynamic Expression Parsing: Evaluates user-defined mathematical strings (e.g., x^3 - 2*x - 5) and computes natural derivatives symbolically on the fly.
2) Step-by-Step State Control: A transport controller (play, pause, next, previous) allows users to scrub through algorithm iterations.
3) Edge-Case Detection: Actively catches and visualizes mathematical breakdowns, such as zero-derivative failures in the Newton-Raphson method or identical interpolation points.
4) Convergence Analytics: Generates a real-time data grid displaying current estimates, function values, and relative error percentages per step.
5) Fully Client-Side: No backend dependencies or database requirements. Computations run locally in the browser.

**Supported Algorithms:**
1) Newton-Raphson Method: Visualizes root-finding via tangent line intersections using symbolic derivatives.
2) Secant Method: Demonstrates root approximation using secant lines between two initial guesses.
3) Muller's Method: Fits a parabola through three points to locate real and complex roots.
4) Natural Cubic Spline Interpolation: Constructs a smooth, piecewise cubic polynomial through a sorted set of user-defined data points.

**Technology Stack:**
1) Frontend: Vanilla HTML5, CSS3, JavaScript (ES6)
2) Canvas Rendering: Native HTML5 <canvas> API (No heavy charting libraries)
3) Mathematics: Math.js (Utilized strictly for safe string parsing and symbolic differentiation)

**How to use the visualizer:**
1) Clone the repository:
2) Navigate to the project directory:
3) Open index.html directly in any modern web browser.

**Author:**
Mantavya Bhojani
