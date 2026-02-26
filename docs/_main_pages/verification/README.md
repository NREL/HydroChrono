# Verification section assets

Verification case images live in one folder per case, so the structure is easy to extend:

- `images_sphere/` — Sphere Model (IEA OES Task 10) verification images
- `images_oswec/` — OSWEC verification images
- `images_rm3_mooring/` — RM3 Mooring (WEC-Sim/MoorDyn) verification images

Pages reference them as `{{ site.baseurl }}/verification/images_<case>/<filename>`.

To add a new verification case: add a new `images_<case_name>/` directory here and place images in it; then link from the case’s verification page using the same URL pattern.
