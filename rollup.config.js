import resolve from '@rollup/plugin-node-resolve';
import { terser } from 'rollup-plugin-terser';

export default {
  input: 'src/beachball.js',
  output: [
    {
      file: 'dist/beachball.js',
      format: 'umd',
      name: 'Beachball'
    },
    {
      file: 'dist/beachball.min.js',
      format: 'umd',
      name: 'Beachball',
      plugins: [terser()]
    }
  ],
  plugins: [resolve()]
};
