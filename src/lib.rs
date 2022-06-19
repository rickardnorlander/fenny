//! # Fenny
//!
//! Fenny is a library for Fenwick Trees, also known as Binary Indexed
//! Trees (BIT). They are a data structure that support three
//! operations in O(log n) time. The first two correspond to
//!
//! ```rust
//! # let (mut a, i, v) = ([0, 0], 1, 1);
//! a[i] += v;
//! # let b:i32 =
//! a[..=i].iter().sum();
//! ```
//!
//! The functions for doing these are respectively `update`,
//! `psum`. There is also `first_larger`, which does a binary search
//! to find the _first_ index for which the prefix sum is _larger_
//! than a some provided value. Being a binary search, it only works
//! if all `a[i]>=0`.
//!
//! ```rust
//! // An array of zeros is a Fenwick tree of zeros.
//! let mut tree = [0; 13];
//! let root = fenny::get_root_p1(13);
//! fenny::update(&mut tree, 2, 7);
//! fenny::update(&mut tree, 5, 20);
//! assert_eq!(fenny::psum(&tree, 2), 7);
//! assert_eq!(fenny::first_larger(&tree, root, 27-1), Some(5));
//! ```
//!
//! ## Flavors
//!
//! The aforementioned functions come in a whole lot of flavors. In
//! addition to the plain ones there are _2d, _3d, _so, _so_2d_lex and
//! _so_3d_lex. You must query using the same flavor that you update
//! with.
//!
//! ### Slope-offset
//!
//! > ⚠️ Slope-offset functions use negative numbers and
//! multiplications internally, select your value-type with this in
//! mind.
//!
//! The _so flavor uses two Fenwick trees. That lets it
//! support the _range update_.
//!
//! ```no_compile
//! a[i..=j] += v;
//! ```
//!
//! which is performed like so
//!
//! ```rust
//! # let (size, i, j, v) = (13, 0, 1, 1);
//! let mut slope = vec![0; size];
//! let mut offset = vec![0; size];
//! fenny::update_so(&mut slope, &mut offset, i..=j, v);
//! ```
//!
//! ### 2d / 3d
//!
//! A 2d (3d) tree is used to represent a 2d(3d) array. The prefix sum
//! for a given point p is over all points in the box defined by p and
//! the origin. These flavors don't support first_larger.  The runtime
//! of update and psum are
//! <nobr>log(size_z) * log(size_y) * log(size_x).</nobr>
//!
//! ```rust
//! use fenny::*;
//! let dim = Dim2{y: 5, x: 8};
//! let mut tree = vec![0; dim.y * dim.x];
//! update_2d(&mut tree, dim , Point2{y: 3, x: 6}, 7);
//! assert_eq!(psum_2d(&tree, dim, Point2{y: 2, x: 6}), 0);
//! assert_eq!(psum_2d(&tree, dim, Point2{y: 3, x: 5}), 0);
//! assert_eq!(psum_2d(&tree, dim, Point2{y: 3, x: 6}), 7);
//! assert_eq!(psum_2d(&tree, dim, Point2{y: 4, x: 7}), 7);
//! ```
//!
//! ### Slope-offset, 2d/3d, lexicographic
//!
//! These flavor uses one tree per dimension to store slope, and one
//! additional tree to store the offset. The uppdate functions operate
//! on a box defined by two points. The psum on the other hand sums in
//! lexicographic order. As first_larger looks for a large enough
//! psum, it too operates in lexicographic order. The runtime of
//! update, psum and first_larger are
//! <nobr>log(size_z) * log(size_y) * log(size_x).</nobr>
//!
//! ```rust
//! use fenny::*;
//! let dim = Dim3{z: 12, y: 6, x: 8};
//! let mut slope_z = vec![0i64; dim.z];
//! let mut slope_y = vec![0i64; dim.z * dim.y];
//! let mut slope_x = vec![0i64; dim.z * dim.y * dim.x];
//! let mut offset  = vec![0i64; dim.z * dim.y * dim.x];
//! let p0 = Point3{z: 3, y: 2, x: 1};
//! let p1 = Point3{z: 7, y: 5, x: 5};
//! update_so_3d_lex(&mut slope_z, &mut slope_y, &mut slope_x, &mut offset, dim, p0, p1, 11);
//! // (4, 0, 0) comes after (3, 2, 1), (3, 2, 2), ... (3, 5, 5)
//! let p2 = Point3{z: 4, y: 0, x: 0};
//! assert_eq!(psum_so_3d_lex(&slope_z, &slope_y, &slope_x, &offset, dim, p2), 5 * 4 * 11);
//! // (4, 0, 0) has a psum of 5*4*11, but so does (3, 5, 5).
//! assert_eq!(first_larger_so_3d_lex(&slope_z, &slope_y, &slope_x, &offset, dim,
//!                                   5 * 4 * 11 - 1), Some(Point3{z:3, y: 5, x: 5}));
//! ```
//!
mod fenny;
pub use crate::fenny::*;
