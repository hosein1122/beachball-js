// __tests__/xy2patch.test.js
import { xy2patch } from "../src/math/geometry.js";

describe("xy2patch", () => {
    test("scalar res scales both axes equally", () => {
        const x = [0, 1, 2];
        const y = [0, 10, 20];
        const res = 2;
        const xy = [5, 5];

        const { vertices, res: resOut, center } = xy2patch(x, y, res, xy);

        expect(vertices).toEqual([
            [5, 5],          // 0*2+5,  0*2+5
            [7, 25],         // 1*2+5, 10*2+5
            [9, 45],         // 2*2+5, 20*2+5
        ]);
        expect(resOut).toEqual([2, 2]);
        expect(center).toEqual([5, 5]);
    });

    test("array res scales X and Y separately", () => {
        const x = [1, 2];
        const y = [3, 4];
        const res = [10, 20];
        const xy = [100, 200];

        const { vertices, res: resOut, center } = xy2patch(x, y, res, xy);

        expect(vertices).toEqual([
            [1 * 10 + 100, 3 * 20 + 200],
            [2 * 10 + 100, 4 * 20 + 200],
        ]);
        expect(resOut).toEqual([10, 20]);
        expect(center).toEqual([100, 200]);
    });

    test("invalid res length throws Error", () => {
        const x = [0];
        const y = [0];
        expect(() => {
            xy2patch(x, y, [2], [0, 0]);
        }).toThrow("res must contain exactly two elements");
    });
});
