import { NodalPlane } from '../src/classes.js';
import { plotDC } from '../src/beachball.js';

describe('plotDC', () => {
    test('returns two patches with vertices', () => {
        const plane = new NodalPlane(30, 45, 90);
        const { tension, pressure } = plotDC(plane, { radius: 100, center: [0, 0] });

        // Basic sanity checks
        expect(tension.vertices.length).toBeGreaterThan(10);
        expect(pressure.vertices.length).toBeGreaterThan(10);

        // Each vertex must be [x,y]
        tension.vertices.forEach(v => expect(v.length).toBe(2));
        pressure.vertices.forEach(v => expect(v.length).toBe(2));

        // Codes must align with vertices array length
        expect(tension.codes.length).toBe(tension.vertices.length);
        expect(pressure.codes.length).toBe(pressure.vertices.length);
    });
});
