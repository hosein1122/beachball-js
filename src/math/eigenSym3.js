/**
 * eigenSym3(A, order = 'asc')
 * @param {number[][]} A    3×3 symmetric matrix
 * @param {'asc' | 'desc'} order  Sort order for eigenvalues and eigenvectors
 * @returns {{values: number[], vectors: number[][]}} Sorted eigenvalues and corresponding eigenvectors
 */
export function eigenSym3(A, order = 'asc') {
    // A is a symmetric 3×3 matrix: [[a11,a12,a13],[a12,a22,a23],[a13,a23,a33]]

    /**
     * Normalize vector sign so the element of largest absolute value is non-negative
     */
    function normalizeSign(v) {
        // find index of max absolute component
        const maxIdx = v.reduce((iMax, x, i, arr) =>
            Math.abs(x) > Math.abs(arr[iMax]) ? i : iMax, 0);
        return v[maxIdx] < 0 ? v.map(x => -x) : v;
    }

    /**
     * Sort eigenvalues and eigenvectors according to the specified order
     */
    function sort(values, vectors) {
        // For ascending: [λ_min, λ_mid, λ_max] => idx = [2,1,0]
        // For descending: [λ_max, λ_mid, λ_min] => idx = [0,1,2]
        const idx = order === 'asc' ? [2, 1, 0] : [0, 1, 2];
        const sortedValues = idx.map(i => values[i]);
        const sortedVectors = idx.map(i => normalizeSign(vectors[i]));
        return { values: sortedValues, vectors: sortedVectors };
    }

    // 1) Scale matrix to avoid overflow
    let s = 0;
    for (let i = 0; i < 3; ++i) {
        for (let j = i; j < 3; ++j) {
            s = Math.max(s, Math.abs(A[i][j]));
        }
    }
    if (s === 0) {
        // Zero matrix -> zero eigenvalues and standard basis eigenvectors
        return sort([0, 0, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    }
    const M = A.map(row => row.map(val => val / s));

    // 2) Compute eigenvalues using the analytical method (Kopp)
    const [m11, m12, m13] = [M[0][0], M[0][1], M[0][2]];
    const [m22, m23, m33] = [M[1][1], M[1][2], M[2][2]];
    const p1 = m12 * m12 + m13 * m13 + m23 * m23;

    if (p1 === 0) {
        // Diagonal matrix case
        const diagVals = [m11, m22, m33].map(v => v * s);
        return sort(diagVals, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    }
    const q = (m11 + m22 + m33) / 3;
    const p2 = (m11 - q) ** 2 + (m22 - q) ** 2 + (m33 - q) ** 2 + 2 * p1;
    const p = Math.sqrt(p2 / 6);

    // Compute B = (M - q I) / p
    const B11 = (m11 - q) / p, B12 = m12 / p, B13 = m13 / p;
    const B22 = (m22 - q) / p, B23 = m23 / p, B33 = (m33 - q) / p;

    // Compute the determinant-like quantity r
    const r = (
        B11 * B22 * B33 + 2 * B12 * B13 * B23
        - B11 * B23 * B23 - B22 * B13 * B13 - B33 * B12 * B12
    ) / 2;
    const rClamped = Math.max(-1, Math.min(1, r));
    const phi = Math.acos(rClamped) / 3;

    // Eigenvalues of the normalized matrix
    const eig1 = q + 2 * p * Math.cos(phi);               // largest
    const eig3 = q + 2 * p * Math.cos(phi + 2 * Math.PI / 3); // smallest
    const eig2 = 3 * q - eig1 - eig3;                     // mid
    const values = [eig1, eig2, eig3].map(v => v * s);

    // 3) Compute eigenvectors via cross-product method
    function eigenvectorFor(lambda) {
        const m = [
            [M[0][0] - lambda, M[0][1], M[0][2]],
            [M[1][0], M[1][1] - lambda, M[1][2]],
            [M[2][0], M[2][1], M[2][2] - lambda]
        ];
        let v = [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0]
        ];
        let norm = Math.hypot(...v);
        if (norm === 0) {
            v = [
                m[0][1] * m[2][2] - m[0][2] * m[2][1],
                m[0][2] * m[2][0] - m[0][0] * m[2][2],
                m[0][0] * m[2][1] - m[0][1] * m[2][0]
            ];
            norm = Math.hypot(...v);
            if (norm === 0) { v = [0, 0, 0]; norm = 1; }
        }
        return v.map(x => x / norm);
    }

    const rawVectors = values.map(val => eigenvectorFor(val / s));

    // 4) Return sorted (and LAPACK‐style normalized) eigenpairs
    return sort(values, rawVectors);
}
