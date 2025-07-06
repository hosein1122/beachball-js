
// beach.js :

import { MomentTensor, NodalPlane } from "../classes";
import { mt2axes } from "./mt2axes";
import { mt2plane } from "./mt2plane";
import { plotDC } from "./plotDC";
import { plotMT } from "./plotMT";


/**
 * Build a generic, renderer-agnostic “collection” for a beach-ball plot.
 * No direct drawing here—just data and style metadata.
 *
 * @param {Array|MomentTensor|NodalPlane} fm
 *   Focal mechanism (3- or 6-element array, or MomentTensor/NodalPlane instance).
 * @param {Object} [opts]
 * @param {number} [opts.linewidth=2]
 * @param {string} [opts.facecolor='b']
 * @param {string} [opts.bgcolor='w']
 * @param {string} [opts.edgecolor='k']
 * @param {number} [opts.alpha=1.0]
 * @param {Array<number>} [opts.xy=[0,0]]
 * @param {number|Array<number>} [opts.width=200]
 * @param {number} [opts.size=100]
 * @param {boolean} [opts.nofill=false]
 * @param {number} [opts.zorder=100]
 * @param {*} [opts.axes=null]   // kept for API compatibility, unused here
 *
 * @returns {Object}
 *   • patches: Array of data-only patch descriptors (ellipse or polygon)  
 *   • fillColors: Array of fill colors (or null if nofill)  
 *   • edgecolor, linewidth, alpha, zorder  
 *   • xy, width, size  
 *   • plotDcUsed: boolean flag which algorithm was used
 */
export function beach(fm, {
    // linewidth = 2,//--
    facecolor = 'b',
    bgcolor = 'w',
    // edgecolor = 'k',//--
    // alpha = 1.0,//--
    xy = [0, 0],
    width = 200,
    size = 100,
    nofill = false,
    // zorder = 100//--
} = {}) {
    // ─── 1) Ensure minimum resolution ───
    if (size < 100) size = 100;

    // ─── 2) Normalize width to [wx, wy] ───
    if (!Array.isArray(width) || width.length !== 2) {
        width = [width, width];
    }

    // ─── 3) Parse fm into np1 and optional MomentTensor mt ───
    let mt = null, np1 = null;
    if (fm instanceof MomentTensor) {
        mt = fm;
        np1 = mt2plane(mt); 
    } else if (fm instanceof NodalPlane) {
        np1 = fm;
    } else if (Array.isArray(fm) && fm.length === 6) {
        mt = new MomentTensor(...fm, 0);
        np1 = mt2plane(mt);
    } else if (Array.isArray(fm) && fm.length === 3) {
        np1 = new NodalPlane(...fm);
    } else {
        throw new TypeError("Wrong input value for 'fm'.");
    }

    // ─── 4) Choose DC or MT algorithm ───
    let colors, patches, plotDcUsed = true;
    if (mt) {
        const [t, n, p] = mt2axes(mt.normalized); 
        if (Math.abs(n.val) < EPSILON && Math.abs(t.val + p.val) < EPSILON) {
            [colors, patches] = plotDC(np1, size, xy, width);
        } else {
            [colors, patches] = plotMT(t, n, p, size,
                true, 0, 0, xy, width);

            plotDcUsed = false;
        }
    } else {
        [colors, patches] = plotDC(np1, size, xy, width);
    }



    // ─── 5) Build fillColors array ───
    const fillColors = patches.map((_, i) => {
        if (nofill)
            return null;
        return colors[i] === 'b' ? facecolor : bgcolor;
    });

    // ─── 6) Return data-only collection ───
    return {
        patches,        // [{ type:'ellipse', xy:[x,y], width, height } | { vertices:[[x,y],…], res:[rx,ry], center:[cx,cy] }]
        fillColors,     // [color|string|null]
        // edgecolor,      // string
        // linewidth,      // number
        // alpha,          // number
        // zorder,         // number
        xy,             // [x,y]
        width,          // [wx,wy]
        size,           // number
        plotDcUsed      // boolean
    };
}
