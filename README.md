# Refugee Gif Viz

This is a visualization of all the refugee immigration data to the US from 2002-2017 via the Refugee Processing Center. The code uses Svelte and D3, both projects by current and former NYT graphics editors, to map lat/lon to SVG and animate things fluidly without introducing unnecessary dependencies.

See the visualization running here (will skip every two weeks and write a frame of a GIF, outputting the final at the end. WARNING: really slow): https://refugee-viz-v0.netlify.com/

## Get started

Install the dependencies...

```bash
cd refugee_gif_viz
npm install
```

...then start [Rollup](https://rollupjs.org):

```bash
npm run dev
```

Navigate to [localhost:5000](http://localhost:5000). You should see your app running. Edit a component file in `src`, save it, and reload the page to see your changes.
