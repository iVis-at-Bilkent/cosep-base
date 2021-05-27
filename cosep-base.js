(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory(require("cose-base"));
	else if(typeof define === 'function' && define.amd)
		define(["cose-base"], factory);
	else if(typeof exports === 'object')
		exports["cosepBase"] = factory(require("cose-base"));
	else
		root["cosepBase"] = factory(root["coseBase"]);
})(self, function(__WEBPACK_EXTERNAL_MODULE__281__) {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ 45:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



var cosepBase = {};

cosepBase.coseBase = __webpack_require__(281);
cosepBase.CoSEPConstants = __webpack_require__(266);
cosepBase.CoSEPEdge = __webpack_require__(466);
cosepBase.CoSEPGraph = __webpack_require__(657);
cosepBase.CoSEPGraphManager = __webpack_require__(303);
cosepBase.CoSEPLayout = __webpack_require__(527);
cosepBase.CoSEPNode = __webpack_require__(794);
cosepBase.CoSEPPortConstraint = __webpack_require__(655);

module.exports = cosepBase;

/***/ }),

/***/ 266:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var CoSEConstants = __webpack_require__(281).CoSEConstants;

function CoSEPConstants() {}

//CoSEPConstants inherits static props in FDLayoutConstants
for (var prop in CoSEConstants) {
  CoSEPConstants[prop] = CoSEConstants[prop];
}

// Initial cooling factors
CoSEPConstants.PHASE2_INITIAL_COOLING_FACTOR = 0.7;
CoSEPConstants.PHASE3_INITIAL_COOLING_FACTOR = 0.5;

// Max iterations for each phase
CoSEPConstants.PHASE1_MAX_ITERATIONS = 2500;
CoSEPConstants.PHASE2_MAX_ITERATIONS = 2500;
CoSEPConstants.PHASE3_MAX_ITERATIONS = 2500;

// Prevent layout from finishing too early for both phases
CoSEPConstants.NOT_TOO_EARLY = 200;

// Default number of ports on one side of a node
CoSEPConstants.PORTS_PER_SIDE = 5;

// # of iterations to check for edge end shifting
CoSEPConstants.EDGE_END_SHIFTING_PERIOD = 5;

// # of iterations to check for node rotation
CoSEPConstants.NODE_ROTATION_PERIOD = 15;

// Thresholds for Phase II
CoSEPConstants.EDGE_END_SHIFTING_FORCE_THRESHOLD = 1;
CoSEPConstants.NODE_ROTATION_FORCE_THRESHOLD = 20;
CoSEPConstants.ROTATION_180_RATIO_THRESHOLD = 0.5;
CoSEPConstants.ROTATION_180_ANGLE_THRESHOLD = 130;

// Polishing (Phase III) Force Constants
CoSEPConstants.DEFAULT_POLISHING_FORCE_STRENGTH = 0.1;

// Further Handling of 1-Degree Nodes
CoSEPConstants.FURTHER_HANDLING_ONE_DEGREE_NODES = true;
CoSEPConstants.FURTHER_HANDLING_ONE_DEGREE_NODES_PERIOD = 50;

module.exports = CoSEPConstants;

/***/ }),

/***/ 466:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var CoSEEdge = __webpack_require__(281).CoSEEdge;
var IGeometry = __webpack_require__(281).layoutBase.IGeometry;
var IMath = __webpack_require__(281).layoutBase.IMath;
var RectangleD = __webpack_require__(281).layoutBase.RectangleD;

function CoSEPEdge(source, target, vEdge) {
    CoSEEdge.call(this, source, target, vEdge);

    // These hold the port constraint for their endpoints.
    // They are divided for easy access
    this.sourceConstraint = null;
    this.targetConstraint = null;
}

CoSEPEdge.prototype = Object.create(CoSEEdge.prototype);

for (var prop in CoSEEdge) {
    CoSEPEdge[prop] = CoSEEdge[prop];
}

// -----------------------------------------------------------------------------
// Section: Getter
// -----------------------------------------------------------------------------

CoSEPEdge.prototype.getSourceConstraint = function () {
    return this.sourceConstraint;
};

CoSEPEdge.prototype.getTargetConstraint = function () {
    return this.targetConstraint;
};

// -----------------------------------------------------------------------------
// Section: Methods
// -----------------------------------------------------------------------------

/**
 * General flag indicating if this edge has any port constraints
 * @returns {boolean}
 */
CoSEPEdge.prototype.isPortConstrainedEdge = function () {
    return !!(this.sourceConstraint || this.targetConstraint);
};

/**
 * Redirects the call to its ports (if any)
 */
CoSEPEdge.prototype.initialPortConfiguration = function () {
    if (!this.isPortConstrainedEdge()) return;

    if (this.sourceConstraint) this.sourceConstraint.initialPortConfiguration();

    if (this.targetConstraint) this.targetConstraint.initialPortConfiguration();

    if (this.targetConstraint && this.sourceConstraint) {
        this.targetConstraint.otherPortConstraint = this.sourceConstraint;
        this.sourceConstraint.otherPortConstraint = this.targetConstraint;
    }
};

/**
 * Changes the calc of edge length based on ports.
 */
CoSEPEdge.prototype.updateLengthWithPorts = function () {
    // If both ends are port constrained then calculate the euler distance between ports
    if (this.sourceConstraint && this.targetConstraint) {
        // If nodes intersect do nothing
        if (this.target.getRect().intersects(this.source.getRect())) {
            this.isOverlapingSourceAndTarget = true;
            return;
        }

        this.isOverlapingSourceAndTarget = false;

        var portSourcePoint = this.sourceConstraint.portLocation;
        var portTargetPoint = this.targetConstraint.portLocation;
        this.lengthX = portTargetPoint.x - portSourcePoint.x;
        this.lengthY = portTargetPoint.y - portSourcePoint.y;
        if (Math.abs(this.lengthX) < 1.0) this.lengthX = IMath.sign(this.lengthX);
        if (Math.abs(this.lengthY) < 1.0) this.lengthY = IMath.sign(this.lengthY);

        this.length = Math.sqrt(this.lengthX * this.lengthX + this.lengthY * this.lengthY);
    }
    // Otherwise, the edge is between one port to a clipping point
    else {
            var clipPointCoordinates = new Array(4);

            if (this.sourceConstraint) {
                this.isOverlapingSourceAndTarget = IGeometry.getIntersection(this.target.getRect(), new RectangleD(this.sourceConstraint.portLocation.x, this.sourceConstraint.portLocation.y, 0, 0), clipPointCoordinates);
            }
            if (this.targetConstraint) {
                this.isOverlapingSourceAndTarget = IGeometry.getIntersection(new RectangleD(this.targetConstraint.portLocation.x, this.targetConstraint.portLocation.y, 0, 0), this.source.getRect(), clipPointCoordinates);
            }

            if (!this.isOverlapingSourceAndTarget) {
                this.lengthX = clipPointCoordinates[0] - clipPointCoordinates[2];
                this.lengthY = clipPointCoordinates[1] - clipPointCoordinates[3];

                if (Math.abs(this.lengthX) < 1.0) this.lengthX = IMath.sign(this.lengthX);

                if (Math.abs(this.lengthY) < 1.0) this.lengthY = IMath.sign(this.lengthY);

                this.length = Math.sqrt(this.lengthX * this.lengthX + this.lengthY * this.lengthY);
            }
        }
};

/**
 * Redirects the call to its ports (if any).
 * Note: Direction is important
 */
CoSEPEdge.prototype.storeRotationalForce = function (springForceX, springForceY) {
    if (!this.isPortConstrainedEdge()) {
        return;
    }

    if (this.sourceConstraint) this.sourceConstraint.storeRotationalForce(springForceX, springForceY);

    if (this.targetConstraint) this.targetConstraint.storeRotationalForce(-springForceX, -springForceY);
};

module.exports = CoSEPEdge;

/***/ }),

/***/ 657:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var CoSEGraph = __webpack_require__(281).CoSEGraph;

function CoSEPGraph(parent, graphMgr, vGraph) {
    CoSEGraph.call(this, parent, graphMgr, vGraph);
}

CoSEPGraph.prototype = Object.create(CoSEGraph.prototype);

for (var prop in CoSEGraph) {
    CoSEPGraph[prop] = CoSEGraph[prop];
}

module.exports = CoSEPGraph;

/***/ }),

/***/ 303:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var CoSEGraphManager = __webpack_require__(281).CoSEGraphManager;

function CoSEPGraphManager(layout) {
    CoSEGraphManager.call(this, layout);

    // All edges with port constraints in this graph manager. The references are hold for efficiency purposes
    this.edgesWithPorts = [];

    // All nodes incident to a port constrained edge.
    this.nodesWithPorts = [];

    // All port constraint endpoints
    this.portConstraints = [];
}

CoSEPGraphManager.prototype = Object.create(CoSEGraphManager.prototype);

for (var prop in CoSEGraphManager) {
    CoSEPGraphManager[prop] = CoSEGraphManager[prop];
}

/**
 * This function needs to update port locations as well.
 */
CoSEPGraphManager.prototype.updateBounds = function () {
    this.rootGraph.updateBounds(true);

    for (var i = 0; i < this.nodesWithPorts.length; i++) {
        this.nodesWithPorts[i].updatePortLocations();
    }
};

module.exports = CoSEPGraphManager;

/***/ }),

/***/ 527:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var FDLayoutConstants = __webpack_require__(281).layoutBase.FDLayoutConstants;
var CoSELayout = __webpack_require__(281).CoSELayout;
var CoSEPConstants = __webpack_require__(266);
var CoSEPGraphManager = __webpack_require__(303);
var CoSEPGraph = __webpack_require__(657);
var CoSEPNode = __webpack_require__(794);
var CoSEPEdge = __webpack_require__(466);

// Constructor
function CoSEPLayout() {
    CoSELayout.call(this);

    // Hold how many ports are available to one side
    this.portsPerSide = CoSEPConstants.PORTS_PER_SIDE;

    // Current phase of the algorithm
    this.phase = CoSEPLayout.PHASE_CORE;
}

CoSEPLayout.prototype = Object.create(CoSELayout.prototype);

for (var property in CoSELayout) {
    CoSEPLayout[property] = CoSELayout[property];
}

// -----------------------------------------------------------------------------
// Section: Class constants
// -----------------------------------------------------------------------------

CoSEPLayout.PHASE_CORE = 1;
CoSEPLayout.PHASE_SECOND = 2;
CoSEPLayout.PHASE_POLISHING = 3;

// -----------------------------------------------------------------------------
// Section: Class methods related to Graph Manager
// -----------------------------------------------------------------------------
CoSEPLayout.prototype.newGraphManager = function () {
    this.graphManager = new CoSEPGraphManager(this);
    return this.graphManager;
};

CoSEPLayout.prototype.newGraph = function (vGraph) {
    return new CoSEPGraph(null, this.graphManager, vGraph);
};

CoSEPLayout.prototype.newNode = function (vNode) {
    return new CoSEPNode(this.graphManager, vNode);
};

CoSEPLayout.prototype.newEdge = function (vEdge) {
    return new CoSEPEdge(null, null, vEdge);
};

// -----------------------------------------------------------------------------
// Section: Limit Phase I CoSE to Specified maxIterations
// -----------------------------------------------------------------------------
CoSEPLayout.prototype.initSpringEmbedder = function () {
    var s = this.getAllNodes().length;
    if (s > FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) {
        this.coolingFactor = Math.max(FDLayoutConstants.COOLING_ADAPTATION_FACTOR, 1.0 - (s - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) / (FDLayoutConstants.ADAPTATION_UPPER_NODE_LIMIT - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) * (1 - FDLayoutConstants.COOLING_ADAPTATION_FACTOR));
    } else {
        this.coolingFactor = 1.0;
    }
    this.initialCoolingFactor = this.coolingFactor;
    this.maxNodeDisplacement = FDLayoutConstants.MAX_NODE_DISPLACEMENT;

    this.maxIterations = CoSEPConstants.PHASE1_MAX_ITERATIONS;

    // Reassign this attribute by using new constant value
    this.displacementThresholdPerNode = 3.0 * FDLayoutConstants.DEFAULT_EDGE_LENGTH / 100;
    this.totalDisplacementThreshold = this.displacementThresholdPerNode * this.getAllNodes().length;

    this.repulsionRange = this.calcRepulsionRange();
};

// -----------------------------------------------------------------------------
// Section: Other Methods
// -----------------------------------------------------------------------------

/**
 * This method introduces port constraints to associated edges. The original CoSE concept is considering edges as a line
 * segment which goes through source node center to target node center but is clipped with respect to node shapes. Now,
 * we want to make sure this point corresponds to a feasible port.
 */
CoSEPLayout.prototype.initialPortConfiguration = function () {
    for (var i = 0; i < this.graphManager.edgesWithPorts.length; i++) {
        var pEdge = this.graphManager.edgesWithPorts[i];
        pEdge.initialPortConfiguration();
    }
};

/**
 * Initialize or reset variables related to the spring embedder
 */
CoSEPLayout.prototype.secondPhaseInit = function () {
    this.phase = CoSEPLayout.PHASE_SECOND;
    this.totalIterations = 0;
    this.maxIterations = CoSEPConstants.PHASE2_MAX_ITERATIONS;
    //this.maxIterations = Math.max(this.getAllNodes().length * 5, CoSEPConstants.PHASE2_MAX_ITERATIONS);

    // Reset variables for cooling
    this.initialCoolingFactor = CoSEPConstants.PHASE2_INITIAL_COOLING_FACTOR;
    this.coolingCycle = 0;
    this.maxCoolingCycle = this.maxIterations / FDLayoutConstants.CONVERGENCE_CHECK_PERIOD;
    this.finalTemperature = FDLayoutConstants.CONVERGENCE_CHECK_PERIOD / this.maxIterations;
    this.coolingAdjuster = 1;

    // Calc of spring forces have to be changes according to ports and stored for edge shifting and rotation
    this.calcSpringForce = function (edge, idealLength) {
        var sourceNode = edge.getSource();
        var targetNode = edge.getTarget();

        // Update edge length
        if (edge.isPortConstrainedEdge()) edge.updateLengthWithPorts();else edge.updateLength();

        if (edge.isOverlapingSourceAndTarget) return;

        var length = edge.getLength();

        if (length == 0) return;

        // Calculate spring forces
        var springForce = edge.edgeElasticity * (length - idealLength);

        // Project force onto x and y axes
        var springForceX = springForce * (edge.lengthX / length);
        var springForceY = springForce * (edge.lengthY / length);

        // Apply forces on the end nodes
        sourceNode.springForceX += springForceX;
        sourceNode.springForceY += springForceY;
        targetNode.springForceX -= springForceX;
        targetNode.springForceY -= springForceY;

        // Store the forces to be used in edge shifting and rotation
        edge.storeRotationalForce(springForceX, springForceY);
    };
};

/**
 * Initialize or reset variables related to the spring embedder
 */
CoSEPLayout.prototype.polishingPhaseInit = function () {
    this.phase = CoSEPLayout.PHASE_POLISHING;
    this.totalIterations = 0;
    this.maxIterations = CoSEPConstants.PHASE3_MAX_ITERATIONS;
    //this.maxIterations = Math.max(this.getAllNodes().length * 5, CoSEPConstants.PHASE3_MAX_ITERATIONS);

    // Node Rotation Related Variables -- No need for rotations/swaps
    // This is for increasing performance in polishing phase
    for (var i = 0; i < this.graphManager.nodesWithPorts.length; i++) {
        var node = this.graphManager.nodesWithPorts[i];
        node.canBeRotated = false;
        node.canBeSwapped = false;
    }

    // Reset variables for cooling
    this.initialCoolingFactor = CoSEPConstants.PHASE3_INITIAL_COOLING_FACTOR;
    this.coolingCycle = 0;
    this.maxCoolingCycle = this.maxIterations / FDLayoutConstants.CONVERGENCE_CHECK_PERIOD;
    this.finalTemperature = FDLayoutConstants.CONVERGENCE_CHECK_PERIOD / this.maxIterations;
    this.coolingAdjuster = 1;
};

/**
 * Here we override the moveNodes method to its FD format again because CoSEP is written 
 * before the changes on the moveNodes method in cose-level are done.
 * So, let's keep its original format to avoid possible effects of the cose-level change.
 */
CoSEPLayout.prototype.moveNodes = function () {
    var lNodes = this.getAllNodes();
    var node;

    for (var i = 0; i < lNodes.length; i++) {
        node = lNodes[i];
        node.move();
    }
};

/**
 * This method implements a spring embedder used by Phase 2 and 3 (polishing) with
 * potentially different parameters.
 *
 * Instead of overloading important functions (e.g. movenodes) we call another fcn so that core CoSE is not affected
 */
CoSEPLayout.prototype.runSpringEmbedderTick = function () {
    this.totalIterations++;

    if (this.totalIterations % CoSEPConstants.CONVERGENCE_CHECK_PERIOD === 0) {
        // If the system is converged
        // But not too early
        if (this.totalIterations >= CoSEPConstants.NOT_TOO_EARLY && this.isConverged()) {
            return true;
        }

        // Update Cooling Temp
        this.coolingCycle++;
        this.coolingAdjuster = this.coolingCycle / 3;
        this.coolingFactor = Math.max(this.initialCoolingFactor - Math.pow(this.coolingCycle, Math.log(100 * (this.initialCoolingFactor - this.finalTemperature)) / Math.log(this.maxCoolingCycle)) / 100 * this.coolingAdjuster, this.finalTemperature);
    }

    this.totalDisplacement = 0;

    // This updates the bounds of compound nodes along with its' ports
    this.graphManager.updateBounds();

    this.calcSpringForces();
    this.calcRepulsionForces();
    this.calcGravitationalForces();
    if (this.phase === CoSEPLayout.PHASE_POLISHING) {
        this.calcPolishingForces();
    }
    this.moveNodes();

    if (this.phase === CoSEPLayout.PHASE_SECOND) {
        if (this.totalIterations % CoSEPConstants.NODE_ROTATION_PERIOD === 0) {
            this.checkForNodeRotation();
        } else if (this.totalIterations % CoSEPConstants.EDGE_END_SHIFTING_PERIOD === 0) {
            this.checkForEdgeShifting();
        }
    }

    if (CoSEPConstants.FURTHER_HANDLING_ONE_DEGREE_NODES && this.totalIterations % CoSEPConstants.FURTHER_HANDLING_ONE_DEGREE_NODES_PERIOD === 0) {
        this.groupOneDegreeNodesAcrossPorts();
    }

    // If we reached max iterations
    return this.totalIterations >= this.maxIterations;
};

/**
 *  Edge shifting during phase II of the algorithm.
 */
CoSEPLayout.prototype.checkForEdgeShifting = function () {
    for (var i = 0; i < this.graphManager.portConstraints.length; i++) {
        this.graphManager.portConstraints[i].checkForEdgeShifting();
    }
};

/**
 * Node Rotation during phase II of the algorithm.
 */
CoSEPLayout.prototype.checkForNodeRotation = function () {
    for (var i = 0; i < this.graphManager.nodesWithPorts.length; i++) {
        this.graphManager.nodesWithPorts[i].checkForNodeRotation();
    }
};

/**
 * Calc polishing forces for ports during phase III (polishing phase)
 */
CoSEPLayout.prototype.calcPolishingForces = function () {
    for (var i = 0; i < this.graphManager.portConstraints.length; i++) {
        this.graphManager.portConstraints[i].calcPolishingForces();
    }
};

/**
 * For one degree nodes incident to port constrained edges, we put them to their desired location if they are not compound nodes
 * and there is one constraint on the edge.
 */
CoSEPLayout.prototype.groupOneDegreeNodesAcrossPorts = function () {
    for (var i = 0; i < this.graphManager.edgesWithPorts.length; i++) {
        var pEdge = this.graphManager.edgesWithPorts[i];
        if (pEdge.sourceConstraint && pEdge.targetConstraint) continue;

        var portConst = pEdge.sourceConstraint || pEdge.targetConstraint;
        if (portConst.otherNode.getEdges().length === 1) {
            var desiredLocation = portConst.getPointOfDesiredLocation();
            portConst.otherNode.setLocation(desiredLocation.getX(), desiredLocation.getY());
        }
    }
};

module.exports = CoSEPLayout;

/***/ }),

/***/ 794:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 */

var CoSEPConstants = __webpack_require__(266);
var CoSENode = __webpack_require__(281).CoSENode;
var PointD = __webpack_require__(281).layoutBase.PointD;
var IMath = __webpack_require__(281).layoutBase.IMath;

function CoSEPNode(gm, loc, size, vNode) {
    CoSENode.call(this, gm, loc, size, vNode);

    // Hold how many ports are available to one side
    this.portsPerSide = CoSEPConstants.PORTS_PER_SIDE;

    // Indicator for having a port constrained edge which the constraint is associated with this node
    this.hasPortConstrainedEdge = false;

    // If the above is true, then this will hold the particular CoSEPPortConstraint classes
    this.associatedPortConstraints = [];

    // In phase II, we will allow nodes with port constrained edges to rotate
    this.canBeRotated = true;

    // In phase II, we will allow nodes with port constrained edges to swap
    this.canBeSwapped = true;

    // Keep rotation history of the node, values to be added: clockwise, counterclockwise, topbottom, rightleft
    // Edit: node rotation and swap are separated later but we continue to use rotationList jointly for both operation
    this.rotationList = [];

    // If the above remains true it will hold the sum of the rotational force induced to this node.
    // Avg can be manually calculated
    this.rotationalForce = 0;

    // This holds the additional force introduced in polishing phase
    this.polishingForceX = 0;
    this.polishingForceY = 0;
}

CoSEPNode.prototype = Object.create(CoSENode.prototype);

for (var prop in CoSENode) {
    CoSEPNode[prop] = CoSENode[prop];
}

// -----------------------------------------------------------------------------
// Section: Methods
// -----------------------------------------------------------------------------

/**
 * This method returns the given port's side and location
 *
 * @param index
 * @returns {any[]}
 *      0 -> Node side
 *      1 -> PointD of port location
 */
CoSEPNode.prototype.getPortCoordinatesByIndex = function (index) {
    var quotient = Math.floor(index / this.portsPerSide);
    var remainder = Math.floor(index % this.portsPerSide);
    var position = void 0;

    switch (quotient) {
        case 0:
            position = new PointD(this.rect.x + this.rect.width * (remainder + 1) / (this.portsPerSide + 1), this.rect.y);
            break;
        case 1:
            position = new PointD(this.rect.x + this.rect.width, this.rect.y + this.rect.height * (remainder + 1) / (this.portsPerSide + 1));
            break;
        case 2:
            position = new PointD(this.rect.x + this.rect.width * (this.portsPerSide - remainder) / (this.portsPerSide + 1), this.rect.y + this.rect.height);
            break;
        case 3:
            position = new PointD(this.rect.x, this.rect.y + this.rect.height * (this.portsPerSide - remainder) / (this.portsPerSide + 1));
            break;
    }

    return [quotient, position];
};

/**
 * This method returns the 'corner' ports of a given side. If there is only one port per side then that port in the only
 * corner port.
 * @param nodeSide
 * @returns {Map<number, [PointD, number]>}
 */
CoSEPNode.prototype.getCornerPortsOfNodeSide = function (nodeSide) {
    var result = new Map();

    if (this.portsPerSide > 1) {
        // Indexes of the ports
        var firstPortIndex = this.portsPerSide * nodeSide;
        var lastPortIndex = this.portsPerSide * (nodeSide + 1) - 1;

        result.set(firstPortIndex, this.getPortCoordinatesByIndex(firstPortIndex)).set(lastPortIndex, this.getPortCoordinatesByIndex(lastPortIndex));
    } else result.set(nodeSide, this.getPortCoordinatesByIndex(nodeSide));

    return result;
};

/**
 * This methods updates all of the associated port constraints locations around itself
 */
CoSEPNode.prototype.updatePortLocations = function () {
    for (var i = 0; i < this.associatedPortConstraints.length; i++) {
        var portConst = this.associatedPortConstraints[i];
        var temp = this.getPortCoordinatesByIndex(portConst.portIndex);
        portConst.portSide = temp[0];
        portConst.portLocation = temp[1];
    }
};

/**
 * Moving a node needs to move its port constraints as well
 *
 * @override
 * @param dx
 * @param dy
 */
CoSEPNode.prototype.moveBy = function (dx, dy) {
    this.rect.x += dx;
    this.rect.y += dy;

    this.associatedPortConstraints.forEach(function (portConstraint) {
        portConstraint.portLocation.x += dx;
        portConstraint.portLocation.y += dy;
    });
};

/**
 * Rotating the node if rotational force inflicted upon is greater than threshold.
 * Sometimes a 180 degree rotation is needed when the above metric does not detect a needed rotation.
 */
CoSEPNode.prototype.checkForNodeRotation = function () {
    if (!this.canBeRotated && !this.canBeSwapped) return;

    // Exceeds threshold? If not then how about 180 degree check for swap (node should also be marked as canBeSwapped)
    var rotationalForceAvg = this.rotationalForce / CoSEPConstants.NODE_ROTATION_PERIOD / this.associatedPortConstraints.length;
    this.rotationalForce = 0;
    if (Math.abs(rotationalForceAvg) < CoSEPConstants.NODE_ROTATION_FORCE_THRESHOLD && this.canBeSwapped) {
        var topBottomRotation = false;
        var rightLeftRotation = false;
        var topBottomPorts = 0;
        var topBottomObstruceAngles = 0;
        var rightLeftPorts = 0;
        var rightLeftObstruceAngles = 0;

        for (var i = 0; i < this.associatedPortConstraints.length; i++) {
            var portConst = this.associatedPortConstraints[i];

            // Don't include one degree nodes to the heuristic
            if (portConst.otherNode.getEdges().length === 1) continue;

            if (portConst.portSide == portConst.sideDirection['Top'] || portConst.portSide == portConst.sideDirection['Bottom']) {
                topBottomPorts++;
                var avgCorrespondingAngle = portConst.correspondingAngle / CoSEPConstants.NODE_ROTATION_PERIOD;
                portConst.correspondingAngle = 0;
                if (avgCorrespondingAngle > CoSEPConstants.ROTATION_180_ANGLE_THRESHOLD) {
                    topBottomObstruceAngles++;
                }
            } else {
                rightLeftPorts++;
                var _avgCorrespondingAngle = portConst.correspondingAngle / CoSEPConstants.NODE_ROTATION_PERIOD;
                portConst.correspondingAngle = 0;
                if (_avgCorrespondingAngle > CoSEPConstants.ROTATION_180_ANGLE_THRESHOLD) {
                    rightLeftObstruceAngles++;
                }
            }
        }

        if (topBottomObstruceAngles / topBottomPorts >= CoSEPConstants.ROTATION_180_RATIO_THRESHOLD) topBottomRotation = true;

        if (rightLeftObstruceAngles / rightLeftPorts >= CoSEPConstants.ROTATION_180_RATIO_THRESHOLD) rightLeftRotation = true;

        if (topBottomRotation) {
            for (var _i = 0; _i < this.associatedPortConstraints.length; _i++) {
                var _portConst = this.associatedPortConstraints[_i];

                if (_portConst.portSide == _portConst.sideDirection['Top']) {
                    _portConst.portIndex = this.portsPerSide - 1 - _portConst.portIndex;
                    _portConst.portIndex = _portConst.portIndex + this.portsPerSide * 2;
                    var temp = this.getPortCoordinatesByIndex(_portConst.portIndex);
                    _portConst.portSide = temp[0];
                    _portConst.portLocation = temp[1];

                    if (_portConst.portConstraintType === _portConst.constraintType['Sided']) for (var _i2 = 0; _i2 < _portConst.portConstraintParameter.length; _i2++) {
                        if (_portConst.portConstraintParameter[_i2] == 0) _portConst.portConstraintParameter[_i2] = 2;
                    }
                } else if (_portConst.portSide == _portConst.sideDirection['Bottom']) {
                    _portConst.portIndex = _portConst.portIndex % this.portsPerSide;
                    _portConst.portIndex = this.portsPerSide - 1 - _portConst.portIndex;
                    var _temp = this.getPortCoordinatesByIndex(_portConst.portIndex);
                    _portConst.portSide = _temp[0];
                    _portConst.portLocation = _temp[1];

                    if (_portConst.portConstraintType === _portConst.constraintType['Sided']) for (var _i3 = 0; _i3 < _portConst.portConstraintParameter.length; _i3++) {
                        if (_portConst.portConstraintParameter[_i3] == 2) _portConst.portConstraintParameter[_i3] = 0;
                    }
                }
            }
            this.rotationList.push("topbottom");
        } else if (rightLeftRotation) {
            for (var _i4 = 0; _i4 < this.associatedPortConstraints.length; _i4++) {
                var _portConst2 = this.associatedPortConstraints[_i4];

                if (_portConst2.portSide == _portConst2.sideDirection['Right']) {
                    _portConst2.portIndex = _portConst2.portIndex % this.portsPerSide;
                    _portConst2.portIndex = this.portsPerSide - 1 - _portConst2.portIndex;
                    _portConst2.portIndex = _portConst2.portIndex + this.portsPerSide * 3;
                    var _temp2 = this.getPortCoordinatesByIndex(_portConst2.portIndex);
                    _portConst2.portSide = _temp2[0];
                    _portConst2.portLocation = _temp2[1];

                    if (_portConst2.portConstraintType === _portConst2.constraintType['Sided']) for (var _i5 = 0; _i5 < _portConst2.portConstraintParameter.length; _i5++) {
                        if (_portConst2.portConstraintParameter[_i5] == 1) _portConst2.portConstraintParameter[_i5] = 3;
                    }
                } else if (_portConst2.portSide == _portConst2.sideDirection['Left']) {
                    _portConst2.portIndex = _portConst2.portIndex % this.portsPerSide;
                    _portConst2.portIndex = this.portsPerSide - 1 - _portConst2.portIndex;
                    _portConst2.portIndex = _portConst2.portIndex + this.portsPerSide;
                    var _temp3 = this.getPortCoordinatesByIndex(_portConst2.portIndex);
                    _portConst2.portSide = _temp3[0];
                    _portConst2.portLocation = _temp3[1];

                    if (_portConst2.portConstraintType === _portConst2.constraintType['Sided']) for (var _i6 = 0; _i6 < _portConst2.portConstraintParameter.length; _i6++) {
                        if (_portConst2.portConstraintParameter[_i6] == 3) _portConst2.portConstraintParameter[_i6] = 1;
                    }
                }
            }
            this.rotationList.push("rightleft");
        }

        return;
    }

    // Exceeds threshold? If yes then rotate node (node should also be marked as canBeRotated)
    if (Math.abs(rotationalForceAvg) >= CoSEPConstants.NODE_ROTATION_FORCE_THRESHOLD && this.canBeRotated) {
        // If this is to be rotated clockwise or counter-clockwise
        var clockwise = Math.sign(rotationalForceAvg) == 1;

        // Change dimension of the node
        var width = this.getWidth();
        var height = this.getHeight();
        this.setWidth(height);
        this.setHeight(width);

        if (clockwise) {
            this.rotationList.push("clockwise");
        } else {
            this.rotationList.push("counterclockwise");
        }

        // Change port locations
        for (var _i7 = 0; _i7 < this.associatedPortConstraints.length; _i7++) {
            var _portConst3 = this.associatedPortConstraints[_i7];

            if (clockwise) {
                _portConst3.portIndex = (_portConst3.portIndex + this.portsPerSide) % (4 * this.portsPerSide);

                if (_portConst3.portConstraintType === _portConst3.constraintType['Sided']) for (var _i8 = 0; _i8 < _portConst3.portConstraintParameter.length; _i8++) {
                    _portConst3.portConstraintParameter[_i8] = (_portConst3.portConstraintParameter[_i8] + 1) % 4;
                }
            } else {
                _portConst3.portIndex = _portConst3.portIndex - this.portsPerSide;
                if (_portConst3.portIndex < 0) _portConst3.portIndex = 4 * this.portsPerSide + _portConst3.portIndex;

                if (_portConst3.portConstraintType === _portConst3.constraintType['Sided']) for (var _i9 = 0; _i9 < _portConst3.portConstraintParameter.length; _i9++) {
                    if (--_portConst3.portConstraintParameter[_i9] < 0) _portConst3.portConstraintParameter[_i9] = 3;
                }
            }

            var _temp4 = this.getPortCoordinatesByIndex(_portConst3.portIndex);
            _portConst3.portSide = _temp4[0];
            _portConst3.portLocation = _temp4[1];
        }
    }
};

/**
 * Modified version of cose to include polishing force
 */
CoSEPNode.prototype.move = function () {
    var layout = this.graphManager.getLayout();
    this.displacementX = layout.coolingFactor * (this.springForceX + this.repulsionForceX + this.gravitationForceX + this.polishingForceX) / this.noOfChildren;

    this.displacementY = layout.coolingFactor * (this.springForceY + this.repulsionForceY + this.gravitationForceY + this.polishingForceY) / this.noOfChildren;

    if (Math.abs(this.displacementX) > layout.coolingFactor * layout.maxNodeDisplacement) {
        this.displacementX = layout.coolingFactor * layout.maxNodeDisplacement * IMath.sign(this.displacementX);
    }

    if (Math.abs(this.displacementY) > layout.coolingFactor * layout.maxNodeDisplacement) {
        this.displacementY = layout.coolingFactor * layout.maxNodeDisplacement * IMath.sign(this.displacementY);
    }

    // a simple node, just move it
    if (this.child == null) {
        this.moveBy(this.displacementX, this.displacementY);
    }
    // an empty compound node, again just move it
    else if (this.child.getNodes().length == 0) {
            this.moveBy(this.displacementX, this.displacementY);
        }
        // non-empty compound node, propogate movement to children as well
        else {
                this.propogateDisplacementToChildren(this.displacementX, this.displacementY);
            }

    layout.totalDisplacement += Math.abs(this.displacementX) + Math.abs(this.displacementY);

    this.springForceX = 0;
    this.springForceY = 0;
    this.repulsionForceX = 0;
    this.repulsionForceY = 0;
    this.gravitationForceX = 0;
    this.gravitationForceY = 0;
    this.polishingForceX = 0;
    this.polishingForceY = 0;
    this.displacementX = 0;
    this.displacementY = 0;
};

/**
 * Overriden here to its old format again 
 * to avoid possible effects of the changes done in cose-level
 * 
 * @param dX
 * @param dY
 */
CoSEPNode.prototype.propogateDisplacementToChildren = function (dX, dY) {
    var nodes = this.getChild().getNodes();
    var node;
    for (var i = 0; i < nodes.length; i++) {
        node = nodes[i];
        if (node.getChild() == null) {
            node.moveBy(dX, dY);
            node.displacementX += dX;
            node.displacementY += dY;
        } else {
            node.propogateDisplacementToChildren(dX, dY);
        }
    }
};

module.exports = CoSEPNode;

/***/ }),

/***/ 655:
/***/ ((module, __unused_webpack_exports, __webpack_require__) => {



/**
 *
 * This object represents the port constraint related to corresponding endpoint
 *
 * @author Alihan Okka
 *
 * @copyright i-Vis Research Group, Bilkent University, 2007 - present
 *
 */

var PointD = __webpack_require__(281).layoutBase.PointD;
var CoSEPConstants = __webpack_require__(266);

function CoSEPPortConstraint(edge, node) {
    // Associated CoSEP Edge
    this.edge = edge;

    // Incident Node of the associated edge
    this.node = node;

    // Holds the type of the constraint using 'enum' object constraintType
    this.portConstraintType = null;

    // Holds additional information related to port's constraint. For instance,
    // Free -> [0, 1, 2, 3]
    // Absolute -> 5    // any numerical value
    // Sided -> [0, 1]   // array of feasible directions
    this.portConstraintParameter = null;

    // Holds the current index at which the port is located at
    this.portIndex = null;

    // Holds the current coordinates of the port
    // { x: - , y: -}
    this.portLocation = null;

    // Holds the direction of the side the port is on
    this.portSide = null;

    // Holds the sum of the rotational force induced to incident node
    // Avg can be manually calculated
    this.rotationalForce = 0;

    // Holds this edges other port constraint (if any)
    this.otherPortConstraint = null;

    // Hold the other node
    this.otherNode = this.edge.getOtherEnd(this.node);

    // Holds the sum of the angle wrt to incident node and desired location.
    // Avg can be manually calculated
    this.correspondingAngle = 0;

    // Hold how many ports are available to one side
    this.portsPerSide = CoSEPConstants.PORTS_PER_SIDE;
}

CoSEPPortConstraint.prototype = Object.create(null);

// -----------------------------------------------------------------------------
// Section: Enumerations
// In some cases, enums are decided implicitly (Top being 0 etc). Thus, be careful when modifying these values
// -----------------------------------------------------------------------------

// An enum to differentiate between different types of constraints
CoSEPPortConstraint.prototype.constraintType = Object.freeze({
    Free: 0,
    Sided: 1,
    Absolute: 2
});

// An enum to  differentiate between different node directions
CoSEPPortConstraint.prototype.sideDirection = Object.freeze({
    Top: 0,
    Right: 1,
    Bottom: 2,
    Left: 3
});

// -----------------------------------------------------------------------------
// Section: Index Iterators
// -----------------------------------------------------------------------------

CoSEPPortConstraint.prototype.nextAdjacentIndex = function () {
    return (this.portIndex + 1) % (4 * this.portsPerSide);
};

CoSEPPortConstraint.prototype.prevAdjacentIndex = function () {
    var temp = this.portIndex - 1;
    if (temp < 0) temp = 4 * this.portsPerSide - 1;

    return temp;
};

CoSEPPortConstraint.prototype.nextAcrossSideIndex = function () {
    return (this.portIndex + 1 + this.portsPerSide) % (4 * this.portsPerSide);
};

CoSEPPortConstraint.prototype.prevAcrossSideIndex = function () {
    var temp = this.portIndex - 1 - this.portsPerSide;
    if (temp < 0) temp = 4 * this.portsPerSide + temp;

    return temp;
};

// -----------------------------------------------------------------------------
// Section: Helper Functions
// -----------------------------------------------------------------------------

// Create a line function going through two given points. Line function returns true if a point is on that line
function line(x1, y1, x2, y2) {
    var slope = (y2 - y1) / (x2 - x1);
    return function (x, y) {
        return y - y1 > slope * (x - x1);
    };
}

// Checks if the testingPoint is left of the line going through point -> otherPoint
function leftTest(point, otherPoint, testingPoint) {
    var test = (otherPoint.x - point.x) * (testingPoint.y - point.y) - (otherPoint.y - point.y) * (testingPoint.x - point.x);

    return test > 0;
}

// -----------------------------------------------------------------------------
// Section: Methods
// -----------------------------------------------------------------------------

/**
 * This method assigns a feasible port to this edge endpoint as follows. For each feasible node side, find the ports
 * closest to node corners. The port with the shortest distance to the other incident node's center is assigned to this
 * port. Obviously, If the port constraint is Absolute, there is nothing to find.
 *
 * Also, add references to CoSEPNode's.
 */
CoSEPPortConstraint.prototype.initialPortConfiguration = function () {
    if (this.portConstraintType === this.constraintType['Absolute']) {
        this.portIndex = this.portConstraintParameter;

        if (this.portIndex > this.portsPerSide * 4 - 1) throw "Error: An absolute port has higher index number than total number of ports!";

        var temp = this.node.getPortCoordinatesByIndex(this.portIndex);
        this.portSide = temp[0];
        this.portLocation = temp[1];

        this.portIndex = this.portConstraintParameter;
    } else {
        // First get all feasible ports
        var allFeasibleCornerPorts = new Map();
        for (var i = 0; i < this.portConstraintParameter.length; i++) {
            var _temp = this.node.getCornerPortsOfNodeSide(this.portConstraintParameter[i]);

            // Merge all feasible ports into allFeasibleCornerPorts
            var _iteratorNormalCompletion = true;
            var _didIteratorError = false;
            var _iteratorError = undefined;

            try {
                for (var _iterator = _temp[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
                    var prop = _step.value;

                    allFeasibleCornerPorts.set(prop[0], prop[1]);
                }
            } catch (err) {
                _didIteratorError = true;
                _iteratorError = err;
            } finally {
                try {
                    if (!_iteratorNormalCompletion && _iterator.return) {
                        _iterator.return();
                    }
                } finally {
                    if (_didIteratorError) {
                        throw _iteratorError;
                    }
                }
            }
        }

        // Find min short distance between ports and other nodes center and assign the port
        var otherNodeCenter = this.otherNode.getCenter();
        var shortestDistance = Number.MAX_SAFE_INTEGER;
        var _iteratorNormalCompletion2 = true;
        var _didIteratorError2 = false;
        var _iteratorError2 = undefined;

        try {
            for (var _iterator2 = allFeasibleCornerPorts[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
                var entry = _step2.value;

                var distance = Math.hypot(entry[1][1].getX() - otherNodeCenter.getX(), entry[1][1].getY() - otherNodeCenter.getY());

                if (distance < shortestDistance) {
                    shortestDistance = distance;
                    this.portIndex = entry[0];
                    this.portSide = entry[1][0];
                    this.portLocation = entry[1][1];
                }
            }
        } catch (err) {
            _didIteratorError2 = true;
            _iteratorError2 = err;
        } finally {
            try {
                if (!_iteratorNormalCompletion2 && _iterator2.return) {
                    _iterator2.return();
                }
            } finally {
                if (_didIteratorError2) {
                    throw _iteratorError2;
                }
            }
        }
    }

    // Adding references
    this.node.hasPortConstrainedEdge = true;
    this.node.associatedPortConstraints.push(this);
    this.node.graphManager.portConstraints.push(this);
};

/**
 * Returns the relative position of port location to related node's center
 * @returns {PointD}
 */
CoSEPPortConstraint.prototype.getRelativeRatiotoNodeCenter = function () {
    var node = this.node;
    return new PointD((this.portLocation.x - node.getCenter().x) / node.getWidth() * 100, (this.portLocation.y - node.getCenter().y) / node.getHeight() * 100);
};

/**
 * The component of the spring force, vertical component if the port is located at the top or bottom and horizontal
 * component otherwise, is considered to be the rotational force.
 *
 *  * The sign of the force should be positive for clockwise, negative for counter-clockwise
 *
 * @param springForceX
 * @param springForceY
 */
CoSEPPortConstraint.prototype.storeRotationalForce = function (springForceX, springForceY) {
    if (this.portSide == this.sideDirection['Top']) {
        this.rotationalForce += springForceX;
        if (this.node.canBeRotated || this.node.canBeSwapped) {
            this.correspondingAngle += Math.abs(this.calcAngle());
            this.node.rotationalForce += springForceX;
        }
    } else if (this.portSide == this.sideDirection['Bottom']) {
        this.rotationalForce -= springForceX;
        if (this.node.canBeRotated || this.node.canBeSwapped) {
            this.correspondingAngle += Math.abs(this.calcAngle());
            this.node.rotationalForce -= springForceX;
        }
    } else if (this.portSide == this.sideDirection['Right']) {
        this.rotationalForce += springForceY;
        if (this.node.canBeRotated || this.node.canBeSwapped) {
            this.correspondingAngle += Math.abs(this.calcAngle());
            this.node.rotationalForce += springForceY;
        }
    } else {
        this.rotationalForce -= springForceY;
        if (this.node.canBeRotated || this.node.canBeSwapped) {
            this.correspondingAngle += Math.abs(this.calcAngle());
            this.node.rotationalForce -= springForceY;
        }
    }
};

/**
 * If the edge is 'Absolute' constrained then there is nothing to do.
 * Otherwise check if the average rotational force inflicted upon port if above threshold.
 * If it exceeds shift the edge (assuming constraint doesn't limit it)
 * Note that there is an additional requirement if the port is located at the 'corner' of node side
 */
CoSEPPortConstraint.prototype.checkForEdgeShifting = function () {
    if (this.portConstraintType == this.constraintType['Absolute']) return;

    // Exceeds threshold?
    // Get AVG and reset the sum
    var rotationalForceAvg = this.rotationalForce / CoSEPConstants.EDGE_END_SHIFTING_PERIOD;
    this.rotationalForce = 0;
    if (Math.abs(rotationalForceAvg) < CoSEPConstants.EDGE_END_SHIFTING_FORCE_THRESHOLD) return;

    var newIndex = null;
    // If the edge wants to go clockwise or counter-clockwise
    if (Math.sign(rotationalForceAvg) == 1) {
        // Currently on a corner port and wants to change node sides. Then check for additional requirements.
        if (this.portIndex % this.portsPerSide == this.portsPerSide - 1) {
            var nextSide = (this.portSide + 1) % 4;
            if (this.portConstraintParameter.includes(nextSide) && this.additionalRequirementForAdjacentSideChanging(nextSide)) {
                newIndex = this.nextAdjacentIndex();
            } else if (this.portConstraintParameter.includes((nextSide + 1) % 4) && this.additionalRequirementForAcrossSideChanging()) {
                newIndex = this.nextAcrossSideIndex();
            }
        } else {
            newIndex = this.nextAdjacentIndex();
        }
    } else {
        if (this.portIndex % this.node.portsPerSide == 0) {
            var _nextSide = this.portSide;
            if (--_nextSide < 0) _nextSide = 3;
            if (this.portConstraintParameter.includes(_nextSide) && this.additionalRequirementForAdjacentSideChanging(_nextSide)) {
                newIndex = this.prevAdjacentIndex();
            } else {
                if (--_nextSide < 0) _nextSide = 3;
                if (this.portConstraintParameter.includes(_nextSide) && this.additionalRequirementForAcrossSideChanging()) {
                    newIndex = this.prevAcrossSideIndex();
                }
            }
        } else {
            newIndex = this.prevAdjacentIndex();
        }
    }

    if (newIndex) {
        var temp = this.node.getPortCoordinatesByIndex(newIndex);
        this.portIndex = newIndex;
        this.portSide = temp[0];
        this.portLocation = temp[1];
    }
};

/**
 * The node needs to be in the right quadrant. Quadrants are defined by node's corner points.
 * Equalities are facing downwards.
 *
 * For line1: Top-Left to Bottom-Right
 * For line2: Top-Right to Bottom-Left
 *
 * @param nextSide
 */
CoSEPPortConstraint.prototype.additionalRequirementForAdjacentSideChanging = function (nextSide) {
    var nodeRect = this.node.rect;
    var otherNodeRect = this.otherNode.getCenter();

    switch (this.portSide) {
        case 0:
            if (nextSide == 1) {
                var check = line(nodeRect.x + nodeRect.width, nodeRect.y, nodeRect.x, nodeRect.y + nodeRect.height);
                return check(otherNodeRect.x, otherNodeRect.y);
            } else {
                var _check = line(nodeRect.x, nodeRect.y, nodeRect.x + nodeRect.width, nodeRect.y + nodeRect.height);
                return _check(otherNodeRect.x, otherNodeRect.y);
            }
        case 1:
            if (nextSide == 0) {
                var _check2 = line(nodeRect.x + nodeRect.width, nodeRect.y, nodeRect.x, nodeRect.y + nodeRect.height);
                return !_check2(otherNodeRect.x, otherNodeRect.y);
            } else {
                var _check3 = line(nodeRect.x, nodeRect.y, nodeRect.x + nodeRect.width, nodeRect.y + nodeRect.height);
                return _check3(otherNodeRect.x, otherNodeRect.y);
            }
        case 2:
            if (nextSide == 3) {
                var _check4 = line(nodeRect.x + nodeRect.width, nodeRect.y, nodeRect.x, nodeRect.y + nodeRect.height);
                return !_check4(otherNodeRect.x, otherNodeRect.y);
            } else {
                var _check5 = line(nodeRect.x, nodeRect.y, nodeRect.x + nodeRect.width, nodeRect.y + nodeRect.height);
                return !_check5(otherNodeRect.x, otherNodeRect.y);
            }
        case 3:
            if (nextSide == 2) {
                var _check6 = line(nodeRect.x + nodeRect.width, nodeRect.y, nodeRect.x, nodeRect.y + nodeRect.height);
                return _check6(otherNodeRect.x, otherNodeRect.y);
            } else {
                var _check7 = line(nodeRect.x, nodeRect.y, nodeRect.x + nodeRect.width, nodeRect.y + nodeRect.height);
                return !_check7(otherNodeRect.x, otherNodeRect.y);
            }
    }
};

/**
 * The node needs to be in the right quadrant. Quadrants are defined by node's center.
 *
 * For line1: Node Center, horizontal, equality pointing downward
 * For line2: Node Center, vertical, equality pointing right
 *
 */
CoSEPPortConstraint.prototype.additionalRequirementForAcrossSideChanging = function () {
    var nodeRect = this.node.getCenter();
    var otherNodeRect = this.otherNode.getCenter();

    switch (this.portSide) {
        case 0:
            return otherNodeRect.y >= nodeRect.y;
        case 1:
            return otherNodeRect.x <= nodeRect.x;
        case 2:
            return otherNodeRect.y <= nodeRect.y;
        case 3:
            return otherNodeRect.x >= nodeRect.x;
    }
};

/**
 * Returns the desired location of the other node/port. It is one idealLength away from this port.
 *
 * @returns {PointD}
 */
CoSEPPortConstraint.prototype.getPointOfDesiredLocation = function () {
    switch (this.portSide) {
        case 0:
            return new PointD(this.portLocation.x, this.portLocation.y - this.edge.idealLength - this.otherNode.getHeight() / 2);
        case 1:
            return new PointD(this.portLocation.x + this.edge.idealLength + this.otherNode.getWidth() / 2, this.portLocation.y);
        case 2:
            return new PointD(this.portLocation.x, this.portLocation.y + this.edge.idealLength + this.otherNode.getHeight() / 2);
        case 3:
            return new PointD(this.portLocation.x - this.edge.idealLength - this.otherNode.getWidth() / 2, this.portLocation.y);
    }
};

/**
 * Calculates the angle between other node/port, this port and desired location of other node/port
 *
 * @returns {number}
 */
CoSEPPortConstraint.prototype.calcAngle = function () {
    var otherPoint = void 0;
    if (this.otherPortConstraint) otherPoint = this.otherPortConstraint.portLocation;else otherPoint = this.otherNode.getCenter();

    var desired = this.getPointOfDesiredLocation();

    var point1 = new PointD(desired.x - this.portLocation.x, desired.y - this.portLocation.y);
    var point2 = new PointD(otherPoint.x - this.portLocation.x, otherPoint.y - this.portLocation.y);

    if (Math.abs(point1.x) < 0) point1.x = 0.0001;
    if (Math.abs(point1.y) < 0) point1.y = 0.0001;

    var angleValue = (point1.x * point2.x + point1.y * point2.y) / (Math.sqrt(point1.x * point1.x + point1.y * point1.y) * Math.sqrt(point2.x * point2.x + point2.y * point2.y));

    var absAngle = Math.abs(Math.acos(angleValue) * 180 / Math.PI);

    return leftTest(this.portLocation, otherPoint, desired) ? -absAngle : absAngle;
};

/**
 * Calculates the polishing force
 */
CoSEPPortConstraint.prototype.calcPolishingForces = function () {
    var edgeVector = new PointD();
    var polishingForceVector = new PointD();

    var otherPoint = void 0;
    if (this.otherPortConstraint) otherPoint = this.otherPortConstraint.portLocation;else otherPoint = this.otherNode.getCenter();

    var desired = this.getPointOfDesiredLocation();

    // Finding the unit vector of the edge
    if (this.edge.getSource() === this.node) {
        edgeVector.setX(this.edge.lengthX / this.edge.length);
        edgeVector.setY(this.edge.lengthY / this.edge.length);
    } else {
        edgeVector.setX(-this.edge.lengthX / this.edge.length);
        edgeVector.setY(-this.edge.lengthY / this.edge.length);
    }

    // Finding which ortogonal unit vector is the one we want
    if (!leftTest(this.portLocation, otherPoint, desired)) {
        polishingForceVector.setX(edgeVector.getY());
        polishingForceVector.setY(-edgeVector.getX());
    } else {
        polishingForceVector.setX(-edgeVector.getY());
        polishingForceVector.setY(edgeVector.getX());
    }

    var distance = Math.hypot(otherPoint.getX() - desired.getX(), otherPoint.getY() - desired.getY());

    var polishingForceX = CoSEPConstants.DEFAULT_POLISHING_FORCE_STRENGTH * polishingForceVector.getX() * distance / 2;
    var polishingForceY = CoSEPConstants.DEFAULT_POLISHING_FORCE_STRENGTH * polishingForceVector.getY() * distance / 2;

    this.otherNode.polishingForceX += polishingForceX;
    this.otherNode.polishingForceY += polishingForceY;
    this.node.polishingForceX -= polishingForceX;
    this.node.polishingForceY -= polishingForceY;
};

module.exports = CoSEPPortConstraint;

/***/ }),

/***/ 281:
/***/ ((module) => {

module.exports = __WEBPACK_EXTERNAL_MODULE__281__;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __webpack_exports__ = __webpack_require__(45);
/******/ 	
/******/ 	return __webpack_exports__;
/******/ })()
;
});