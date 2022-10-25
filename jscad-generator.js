/**
 * Rounded Cuboid
 * @category Creating Shapes
 * @skillLevel 1
 * @description Demonstrating the roundedCuboid() functiom
 * @tags cube, cuboid, shape, 3d, rounded
 * @authors Rene K. Mueller
 * @licence MIT License
 */

const { cuboid, roundedCuboid, ellipsoid, cylinder, cylinderElliptic, star, geodesicSphere } = require('@jscad/modeling').primitives
const { extrudeLinear } = require('@jscad/modeling').extrusions
const { bezier } = require('@jscad/modeling').curves
const { hull, hullChain } = require('@jscad/modeling').hulls
const { measureAggregateArea, measureAggregateVolume } = require('@jscad/modeling').measurements;
const { translate, translateX, translateY, translateZ, scale, rotate, rotateX, rotateY, rotateZ } = require('@jscad/modeling').transforms

 
// example component tree:
const exampleComponentModel = 
  { type: "root", params: {segments: 128}, children: [
    { type: "SingleHull", params: {length: 20, width: 10, height: 7}, children: [
      { type: "Front", children: [
        { type: "LinearNoseEllipsoid", params: { noseLength: 14, endWidth: 10, endHeight: 7 } }
      ]},
      { type: "Midsection", children: [
        { type: "Midshape", children: [
          // { type: "Boxfish", params: { radius: 2 } }
          { type: "EllipseBody" }  // cylinder is a special case where width = height
          // { type: "Star", params: { points: 5, innerRadius: 7, outerRadiusFactor: 2 } }
        ]}
      ]},
      { type: "Tail", children: [
        // { type: "LinearTailCone", params: { tailLength: 3, endRadius: 2 } }
        { type: "LinearTailEllipsoid", params: { tailLength: 24, endWidth: 24, endHeight: 24 } }
      ]}
    ]}
  ]};

const getParameterDefinitions = () => [
  { name: 'jsonComponentModel', type: 'text', default: JSON.stringify(exampleComponentModel), caption: "JSON hierarchy:" },
  { name: 'cx', type: 'float', default: 2000, caption: "Center X" },
  { name: 'cy', type: 'float', default: 500, caption: "Center Y" },
  { name: 'cz', type: 'float', default: 500, caption: "Center Z" },
  // { name: 'length', type: 'float', default: 12, caption: 'Length:' },
  // { name: 'width', type: 'float', default: 7, caption: 'Width:' },
  // { name: 'height', type: 'float', default: 7, caption: 'Height:' },
  // { name: 'nose', type: 'float', default: 3, caption: 'Nose:' },
  // { name: 'tail', type: 'float', default: 3, caption: 'Tail:' },
  // { name: 'rounded', type: 'choice', caption: 'Round the corners', values: [0, 1], captions: ['No', 'Yes'], default: 1 },
  // { name: 'radius', type: 'float', default: 2, caption: 'Corner Radius:' },
  // { name: 'segments', type: 'int', default: 32, caption: "Segments:"}
]


let wingExampleComponentModel = 
  { type: "root", params: {segments: 128}, children: [
    { type: "SingleHull", params: {length: 18, width: 7, height: 7}, children: [
      { type: "Front", children: [
        { type: "LinearNoseCone", params: { noseLength: 3 } }
      ]},
      { type: "Midsection", children: [
        { type: "Wing", params: { wingSpan: 18, wingDepth: 5 } },
        { type: "Midshape", children: [
          { type: "Boxfish", params: { radius: 2} }
        ]}
      ]},
      { type: "Tail", children: [
        { type: "LinearTailCone", params: { tailLength: 3 } }
      ]}
    ]}
  ]};


let exampleComponentGenerator = (str) => {
  // convert a non-json CFG-style string into a component tree.
  // TODO(zamfi): implement this
}

const grammar = {
  root: {params: ["segments"], children: [ "SingleHull" ] },
  SingleHull: {params: ["length", "width", "height"], children: [ ["Front", "Midsection", "Tail"] ] },
  Front: { children: [ "LinearNoseEllipsoid", "LinearNoseCone" ] },
  Midsection: { children: [ "Midshape" /*, "Wing" */ ] },
  Tail: { children: [ "LinearTailEllipsoid", "LinearTailCone" ] },
  Midshape: { children: [ "Boxfish", "EllipseBody", "Star" ] },
  Boxfish: { params: ["radius"], children: null },
  EllipseBody: { children: null },
  Star: { params: ["points", "innerRadius", "outerRadiusFactor"], children: null },
  LinearNoseEllipsoid: { params: ["noseLength", "endWidth", "endHeight"], children: null },
  LinearNoseCone: { params: ["noseLength"], children: null },
  LinearTailEllipsoid: { params: ["tailLength", "endWidth", "endHeight"], children: null },
  LinearTailCone: { params: ["tailLength", "endRadius"], children: null }
}

const parameterRanges = {
  root: { segments: { min: 128, max: 128 } },
  SingleHull: { length: { min: 10, max: 40 }, width: { min: 1, max: 20 }, height: { min: 1, max: 20 } },
  Boxfish: { radius: { min: 0, max: 5 } },
  Star: { points: { min: 3, max: 10 }, innerRadius: { min: 1, max: 10 }, outerRadiusFactor: { min: 1, max: 10 } },
  LinearNoseEllipsoid: { noseLength: { min: 1, max: 10 }, endWidth: { min: 1, max: 10 }, endHeight: { min: 1, max: 10 } },
  LinearNoseCone: { noseLength: { min: 1, max: 10 } },
  LinearTailEllipsoid: { tailLength: { min: 1, max: 10 }, endWidth: { min: 1, max: 10 }, endHeight: { min: 1, max: 10 } },
  LinearTailCone: { tailLength: { min: 1, max: 10 }, endRadius: { min: 1, max: 10 } }
}


const generate = (node, contextParams={}) => {
  let nodeType = node.type;
  let nodeParams = node.params || {};
  let params = {...contextParams, ...nodeParams}; // node params shadow context params
  console.log("generating", nodeType, "with params", params);

  switch (nodeType) {
    // start
    case "root":
      return node.children.flatMap(n => generate(n, params));
    
    // branch
    case "SingleHull":
      return [hull(...node.children.flatMap(n => generate(n, params)))];
    case "Front":
    case "Midsection":
    case "Midshape":
    case "Tail":
      return node.children.flatMap(n => generate(n, params));
    
    // leaves
    case "LinearNoseEllipsoid":
      return generateNoseEllipsoid(params);
    // case "NoseGeodesic":
    //   return generateNoseGeodesic(params);
    case "LinearNoseCone":
      return generateNoseCone(params);
    case "LinearTailCone":
      return generateTailCone(params);
    case "LinearTailEllipsoid":
      return generateTailEllipsoid(params);
    case "Boxfish":
      return generateBoxfish(params);
    case "EllipseBody":
      return generateEllipseBody(params);
    case "Star":  // doesn't actually make a star b/c smoothing algo, just an extruded polygon
      return generateStar(params);

    default:
      return [];
  }
}

const generateBoxfish = (params) => {
  if (params.radius > 0) {
    return [ roundedCuboid({ size: [params.length, params.width, params.height], roundRadius: params.radius, segments: params.segments }) ];
  } else {
    return [ cuboid({ size: [params.length, params.width, params.height] }) ];
  }
}

const generateEllipseBody = (params) => {
  return [ rotateY(Math.PI/2,
    cylinderElliptic({ height: params.length, 
      startRadius: [params.height, params.width], 
      endRadius: [params.height, params.width], 
      segments: params.segments }) 
  )];
}

const generateStar = (params) => {
  return [ rotateY(Math.PI/2,
    extrudeLinear({height: params.length}, 
      star({vertices: params.points, innerRadius: params.innerRadius, outerRadius: params.innerRadius * params.outerRadiusFactor})
      )
  )]
}

const generateNoseCone = (params) => {
  return [
    ellipsoid({radius: [params.length/10, params.width/10, params.height/10],
                              center: [params.length*0.5+params.noseLength, 0, 0],
                              segments: params.segments})
  ];
}

const generateNoseEllipsoid = (params) => {
  return [
    ellipsoid({radius: [params.noseLength, params.endWidth, params.endHeight],
                              center: [params.length*0.5+params.noseLength, 0, 0],
                              segments: params.segments})
  ];
}

// Doesn't work for whatever reason
// const generateNoseGeodesic = (params) => {
//   return [translateX((params.length*(0.5+params.noseLength/10)),
//     geodesicSphere({radius: params.geodesicRadius, frequency: params.geodesicFrequency})
//   )];
// }

const generateTailCone = (params) => {
  return [
    translateX(-(params.length*0.5+params.tailLength),
                                  rotateY(Math.PI/2, 
                                          cylinder({radius: params.endRadius, 
                                                    height: params.tailLength,
                                                    segments: params.segments})))
  ];
}

const generateTailEllipsoid = (params) => {
  return [
    ellipsoid({radius: [params.tailLength, params.endWidth, params.endHeight],
                              center: [-params.length*0.5+params.tailLength, 0, 0],
                              segments: params.segments})
  ];
}

 /**
  * Create a rounded cuboid with the supplied parameters
  * @param {Number} params.width - The cuboid's width.
  * @param {Number} params.depth - The cuboid's depth.
  * @param {Number} params.height - The cuboid's height.
  * @param {Number} params.rounded - 1 if the cuboid should be rounded, 0 otherwise.
  * @param {Number} params.radius - The cuboid's corner radius.
  * @returns {geometry}
  */
const main = (params) => {
  let {cx, cy, cz, jsonComponentModel} = params;
  let componentModel = JSON.parse(jsonComponentModel);
  let generatedModel = generate(componentModel);
  console.log("==fairing metadata==")
  console.log(JSON.stringify({
    area: measureAggregateArea(generatedModel),
    volume: measureAggregateVolume(generatedModel)
  }));
  return translate([cx, cy, cz], rotateZ(Math.PI, generatedModel));
}
 
module.exports = { main, getParameterDefinitions, grammar, parameterRanges }
 
 
 
 