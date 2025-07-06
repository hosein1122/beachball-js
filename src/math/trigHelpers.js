//trigHelpers.js:
//trigonometry Helpers

// Angle conversions
export const D2R = Math.PI / 180;  // degrees → radians
export const R2D = 180 / Math.PI;  // radians → degrees

// Small tolerance value
export const EPSILON = 1e-5;



/**
* Notes:
* - If either input is an array, the other is broadcast to the same length.
* - Throws if array lengths mismatch.
 *
 * @param {number | number[]} th angle(s) in radians
 * @param {number | number[]} r radius/radii
 * @returns {[number | number[], number | number[]]} (same shape as inputs)
 */
export function pol2cart(th, r) {
    const isArrayTH = Array.isArray(th);
    const isArrayR = Array.isArray(r);

    // Helper for scalar calculations
    const calc = (theta, radius) => {
        const xVal = radius * Math.cos(theta);
        const yVal = radius * Math.sin(theta);
        return [xVal, yVal];
    };

    // Case 1: both scalars
    if (!isArrayTH && !isArrayR) {
        return calc(th, r);
    }

    // Case 2: broadcast arrays
    const thArr = isArrayTH ? th : Array(r.length).fill(th);
    const rArr = isArrayR ? r : Array(th.length).fill(r);

    if (thArr.length !== rArr.length) {
        throw new Error('pol2cart: th and r must have the same length when arrays');
    }

    const xs = [];
    const ys = [];
    for (let i = 0; i < thArr.length; i += 1) {
        const [xi, yi] = calc(thArr[i], rArr[i]);
        xs.push(xi);
        ys.push(yi);
    }
    return [xs, ys];
}


