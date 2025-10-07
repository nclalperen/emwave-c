# emwave-c

Real-time 2D FDTD EM simulator (C + SDL2). Cross-platform via CMake.

## Key achievements
- Real-time 512×512 @ 60–120 FPS; 1024×768 @ 30–60 FPS; <16 ms UI latency.
- CPML with presets + manual tuning; reflections < -60 dB at target f.
- Instrumentation: dual probes, scope + FFT, Auto-P99/Peak/Hold, paint/edit tools.
- Portable Windows .exe (~1–2 MB), reproducible builds (Linux/Windows).
