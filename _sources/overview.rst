.. _what:

What is Easy Climate
====================================

.. figure:: _static/fig1.jpg
    :scale: 40%
    :align: center

    Photo by Noam Almosnino on Unsplash

.. warning::

    Due to the continuous development of the repository, some content may become outdated. Please refer to the `API <https://easyclimate.readthedocs.io/en/latest/api_index/index.html>`_ to clarify the functionality of this repository.

Hey there, climate wrangler! 😎 Ever stared at a NetCDF file the size of a small country's ego,
wondering if it's plotting against you? Or wrestled with EOF decompositions that feel like they're decomposing your sanity? 😩
Fear not—EasyClimate is here to swoop in like a caffeinated superhero,
turning those multi-hour code marathons into "one-liner" victory laps.
We're talking Python-powered wizardry that crunches terabytes of gridded goodies—from
ERA5 reanalysis to WRF model outputs—while you sip coffee and ponder the mysteries of El Niño. 🎉

Built for the bold (that's you, developer extraordinaire), EasyClimate isn't just a package;
it's an ecosystem of front-end flair and back-end brawn. Prototype like a poet in Python,
scale like a boss with compiled kernels, and maybe even flirt with Rust for that extra concurrency kick.
Let's dive in—shall we? 🚀


Why EasyClimate? (Because Who Has Time for Tedious?)
------------------------------------------------------------------

Tired of glue-code graveyards and visualization voodoo? 😒 We've got your back with:

- One-Line Wonders: Load, slice, analyze, and plot? Done. No more "import everything under the sun" rituals.
- Speed Demon Mode: Backend offloads the heavy math (think finite differences on steroids) so you don't wait for results like it's 1999 dial-up.
- Modular Magic: Pick 'n' mix modules for stats, physics, filters, and more—extensible like your favorite LEGO set, but for atmospheric science.
- Open-Source Shenanigans: BSD-3-Clause licensed, community-fueled, and begging for your pull requests. Join the party! 👨‍👩‍👧‍👦

Unlock insights into ENSO shenanigans, typhoon tantrums, or ocean mixed-layer moods—faster than you can say "climate attribution." 🌍💨


The Squad: Frontend, Backend, and That Rust Enigma
EasyClimate's like a heist crew: Python's the charming leader, Backend's the muscle, and Rust?
The mysterious wildcard whispering promises of memory safety and parallelism. Here's the lineup:

.. tip::

    Frontend's your daily driver (pip install easyclimate), but hook up Backend for that "whoa, it's fast" moment. Dependencies? NumPy, xarray, SciPy. No Dask drama here—we keep it lean, but scalable. 🏗️


Module Mayhem: What's in the Toolbox?
------------------------------------------------------------------

Import once (``import easyclimate as ecl``) and let the good times roll. We've leveled up from the old docs—now with fresher physics,
field-specific flair, and WRF whispers. Here's the squad, post-2025 glow-up:

- Core: The brainiac hub—data nodes, EOF/MCA (now with backend turbo), trends (Mann-Kendall on steroids), variability vibes, and spherical harmonics for those global groove analyses.
- Physics: Geo-gems (Coriolis curls), moisture mysteries (dewpoint divas, mixing ratios), temp tricks (virtual temps, potential θ), convection chaos (stability checks), and energy equations that won't make you cry.
- Filter: Time tamer (Butterworth butter-ups, Lanczos low-passes, wavelet wiggles, EMD decomps) + spatial sorcery (Barnes blending, Gaussian glows) + spectra sleuthing (Redfit reds, anyone?).
- Interp: Regrid rebels (Barnes objective analysis, mesh-to-mesh morphs, hybrid-to-pressure hops)—interpolate like a pro, no interpolation angst.
- Plot: Visual virtuosos—axis acrobats, projection parties (ortho, polar), Taylor diagrams for correlation confetti, quick-draw xarray sketches, and curved quivers that quiver just right. Cartopy + matplotlib = eye candy. 🎨
- Field: Climate classics—teleconnections (PNA parties, AO auras), monsoons (BSISO beats), ocean ops (MLD depths, front finds), typhoons (intensity intel), waves (MJO marches), and heat stress (WBGT warnings).
- WRF: Model whisperer—post-process those Weather Research & Forecasting outputs like a breeze. Native spherical harmonics via PySpharm/WindSpharm for extra spin.

Full API deets? Hit the `docs <https://easyclimate.readthedocs.io/en/latest/api_index/index.html>`_—we've sprinkled in new gems like equatorial wave spectra and moist lapse rates since the last overview dust-up.


Join EasyClimate: Contribute & Crystal Ball
------------------------------------------------------------------

Got a killer kernel or a quirky index? Fork us on GitHub—Backend begs for ARM/GPU optimizations,
Frontend craves custom fields (heatwaves, anyone?). Roadmap? Stabilize those APIs,
Rust-ify the concurrency, and satellite modules for starry-eyed expansions.

Questions? Hit the repos or Zenodo. Let's make climate code cool again—together.
What's your first hack? Drop it in the issues! 🛠️🎊

EasyClimate Related Repositories
------------------------------------------------------------------

The EasyClimate ecosystem consists of several interconnected repositories:

.. figure:: _static/easyclimate_backend_logo_mini.png
    :scale: 20%
    :align: center
    :target: https://github.com/shenyulu/easyclimate-backend

    The backend of easyclimate (mainly Fortran) (https://github.com/shenyulu/easyclimate-backend)

.. figure:: _static/easyclimate_rust_logo_mini.png
    :scale: 20%
    :align: center
    :target: https://github.com/shenyulu/easyclimate-rust

    The Rust backend of easyclimate (https://github.com/shenyulu/easyclimate-rust)

.. figure:: _static/easyclimate_map_logo_mini.png
    :scale: 20%
    :align: center
    :target: https://github.com/shenyulu/easyclimate-map

    Easily obtain map data (https://github.com/shenyulu/easyclimate-map)
