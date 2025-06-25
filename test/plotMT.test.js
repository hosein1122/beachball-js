import { MomentTensor } from '../src/classes.js';
import { mt2axes, plotMT } from '../src/beachball.js';

describe('plotMT early returns', () => {
    test('pure explosion/implosion yields no segments', () => {
        // All eigenvalues equal → pure isotropic
        const mt = new MomentTensor([1, 1, 1, 0, 0, 0], 0);
        const [T, N, P] = mt2axes(mt);
        const { tension, pressure } = plotMT(T, N, P, { radius: 100, center: [0, 0] });

        expect(Array.isArray(tension.vertices)).toBe(true);
        expect(tension.vertices.length).toBe(0);
        expect(Array.isArray(pressure.vertices)).toBe(true);
        expect(pressure.vertices.length).toBe(0);
    });
});

describe('plotMT general case', () => {
    test('simple double-couple produces segments', () => {
        // A simple tensor with non-zero deviatoric part
        const mt = new MomentTensor([1, 1, -2, 0, 0, 0], 0);
        const [T, N, P] = mt2axes(mt);
        const { tension, pressure } = plotMT(T, N, P, { radius: 100, center: [0, 0] });

        // Expect at least one vertex in each patch
        expect(tension.vertices.length).toBeGreaterThan(0);
        expect(pressure.vertices.length).toBeGreaterThan(0);

        // codes length matches vertices length
        expect(tension.codes.length).toBe(tension.vertices.length);
        expect(pressure.codes.length).toBe(pressure.vertices.length);
    });
});
