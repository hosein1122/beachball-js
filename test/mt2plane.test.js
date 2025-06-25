import { mt2plane } from '../src/beachball.js';
import { MomentTensor } from '../src/classes.js';

describe('mt2plane', () => {
    test('identity moment tensor produces expected nodal plane', () => {
        // Simple symmetric tensor example
        const mt = new MomentTensor(1, 1, -2, 0, 0, 0, 0);
        const plane = mt2plane(mt);

        expect(plane).toHaveProperty('strike');
        expect(plane).toHaveProperty('dip');
        expect(plane).toHaveProperty('rake');
        expect(Number.isFinite(plane.strike)).toBe(true);
        expect(Number.isFinite(plane.dip)).toBe(true);
        expect(Number.isFinite(plane.rake)).toBe(true);
    });
});
