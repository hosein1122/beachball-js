import { auxPlane } from '../src/beachball.js';

describe('auxPlane', () => {
    test('plane (30°, 45°, 90°) → (210°, 45°, 90°)', () => {
        const [strike, dip, rake] = auxPlane(30, 45, 90);
        expect(strike).toBeCloseTo(210);
        expect(dip).toBeCloseTo(45);
        expect(rake).toBeCloseTo(90);
    });
});
