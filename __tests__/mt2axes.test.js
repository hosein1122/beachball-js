// __tests__/mt2axes.test.js
import { mt2axes } from "../src/beachball/mt2axes";
import { MomentTensor, PrincipalAxis } from "../src/classes";

describe("mt2axes", () => {
    test("diagonal moment tensor → axes align with coordinate axes", () => {
        // λmin=1, λmid=2, λmax=3; eigenvectors = standard basis
        let mt = new MomentTensor(1, 2, 3, 0, 0, 0, 0);

        const [T, N, P] = mt2axes(mt);

        // 1) eigenvalues → T=3, N=2, P=1
        expect(T.val).toBeCloseTo(3, 6);
        expect(N.val).toBeCloseTo(2, 6);
        expect(P.val).toBeCloseTo(1, 6);

        // 2) for this tensor:
        //    • T‐axis vector = (0,0,1) → plunge = 0°, strike = 270°
        //    • N‐axis vector = (0,1,0) → plunge = 0°, strike = 360°
        //    • P‐axis vector = (1,0,0) → plunge = 90°, strike = 360°
        expect(T.strike).toBeCloseTo(270, 6);
        expect(T.dip).toBeCloseTo(0, 6);

        expect(N.strike).toBeCloseTo(360, 6);
        expect(N.dip).toBeCloseTo(0, 6);

        expect(P.strike).toBeCloseTo(360, 6);
        expect(P.dip).toBeCloseTo(90, 6);

        // also verify instances
        expect(T).toBeInstanceOf(PrincipalAxis);
        expect(N).toBeInstanceOf(PrincipalAxis);
        expect(P).toBeInstanceOf(PrincipalAxis);
    });
});
