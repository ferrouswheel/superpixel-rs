# Superpixels for Rust

This crate contains algorithms for segmenting an image into perceptually similar regions.

Currently it only has one implementation:
- [SLIC (Simple Linear Iterative Clustering)](http://ivrl.epfl.ch/research/superpixels) translated to Rust.

## Example Command-line Usage

```
cargo run --features="build-binary" -- \
    --outline-image=assets/hawaii_outline.jpg assets/hawaii.jpg
```

![Test image](assets/hawaii.jpg?raw=true "Test image")
![Segmented Test image](assets/hawaii_outline.jpg?raw=true "Segmented Test image")

## Library Usage

The `superpixel` executable is a reasonably compact example. This crate expects
images in the format of the `image` crate, specially a `image::DynamicImage::Rgb8`.


```rust
extern crate superpixel;
extern crate image;

use std::vec::Vec;
use std::iter::repeat;
use std::time::Instant;

use image::DynamicImage;

use superpixel::{LABColor,rgb_to_lab_image,do_segmentation_with_num_superpixels,save_image_with_segment_boundaries};

...

let img : image:DynamicImage::Rgb8 = ...; // get your image from somewhere

// SLIC operates on LAB space
let mut lab_vec : Vec<LABColor> = repeat(Default::default()).take((img.width() * img.height()) as usize).collect();
rgb_to_lab_image(&img, &mut lab_vec);

let num_patches = 9000;
let compactness = 20;
let iterations = 10;

// klabels has the pixel labelling to segment id
let mut klabels : Vec<i32> = vec![];
// while 9000 superpixels/patches are requested, the exact number may differ so the actual number is returned
let num_superpixels = do_segmentation_with_num_superpixels(
    &img, &mut lab_vec,
    img.width(), img.height(),
    &mut klabels,
    num_patches, compactness, iterations);

// debug image for showing where segment boundaries are
let outline_fn = "test.png";
let outline_color = image::Rgb([255u8, 255u8, 255u8]);
save_image_with_segment_boundaries(&img, &klabels, &outline_fn, &outline_color);
```

## Windows

When developing on Windows - if you have installed 64bit libraries, and 64bit MSVC rust, you need to do

```
"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x64
```

before running `cargo build`