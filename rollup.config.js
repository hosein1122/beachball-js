import resolve from '@rollup/plugin-node-resolve';
import commonjs from '@rollup/plugin-commonjs';
import terser from '@rollup/plugin-terser';

const terserOptions = {
  compress: { drop_console: true, passes: 2 },
  mangle: { toplevel: true },
  format: { comments: false }
};

export default {
  input: 'src/beachball-component.js',

  /* ماژول‌های Node-only را بیرونی کنیم */
  external: ['canvas', 'fs', 'path', 'stream', 'buffer', 'util'],

  treeshake: { moduleSideEffects: false },

  output: [
    { file: 'dist/beachball.esm.js', format: 'es' },
    { file: 'dist/beachball.esm.min.js', format: 'es', plugins: [terser(terserOptions)] },
    { file: 'dist/beachball.umd.js', format: 'umd', name: 'Beachball' },
    {
      file: 'dist/beachball.umd.min.js', format: 'umd', name: 'Beachball',
      plugins: [terser(terserOptions)]
    },
  ],

  plugins: [
    resolve({ browser: true, preferBuiltins: false }),
    commonjs(),
  ]
};