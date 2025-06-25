import { beach } from '../src/beachball.js';
import { NodalPlane, MomentTensor } from '../src/classes.js';

describe('beach', () => {
    test('nodal plane input yields two patches', () => {
        const bp = beach([30, 45, 90], { radius: 100 });
        expect(bp.tension.vertices.length).toBeGreaterThan(0);
        expect(bp.pressure.codes.length).toBe(bp.pressure.vertices.length);
        expect(bp.tension.facecolor).toBe('b');
    });

    test('moment tensor input yields styled patches', () => {
        const mt = new MomentTensor(1, 1, -2, 0, 0, 0, 0);
        const bp = beach(mt, { radius: 80, facecolor: 'red', bgcolor: 'blue' });
        expect(bp.tension.facecolor).toBe('red');
        expect(bp.pressure.facecolor).toBe('blue');
    });
});
