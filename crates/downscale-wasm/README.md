# downscale-wasm

WebAssembly bindings for downscale-rs — run bias correction in the browser,
no server, ~74 KB wasm.

## Build

```bash
wasm-pack build --target web --out-dir pkg
```

This generates `pkg/` (gitignored): the `.wasm` binary plus the JS/TS glue.

## Demo

```bash
# from this crate's directory, after building pkg/
python3 -m http.server 8000
# open http://localhost:8000/www/index.html
```

`www/index.html` generates a synthetic "observed" climate and a biased
"model", corrects it with empirical quantile mapping, and plots the three
empirical CDFs with the before/after metrics (KS, mean bias, RMSE).

## API

```js
import init, {
  QuantileMapping, QuantileDeltaMapping, rmse, ksStatistic, meanBias,
} from "./pkg/downscale_wasm.js";

await init();
const qm = new QuantileMapping(obs, model, 100, "mult", "midpoint"); // Float64Array
const corrected = qm.apply(model);                                   // Float64Array
qm.free();                                                           // free wasm memory
```

Series are `Float64Array`; `kind` ∈ `{"add","mult"}`, `nodes` ∈
`{"endpoints","midpoint"}`. Errors surface as JS exceptions.
