# Fenny

Fenny is a library for Fenwick Trees, also known as Binary Indexed
Trees (BIT). They are a data structure that can efficiently (as in
`O(log n)`) support three operations.

```rust
a[i] += v;
a[..=i].iter().sum();
```

The third operation is a binary search on prefix sums.

Here's a simple example usage:

```rust
let mut tree = [0; 13];   // An array of zeros is a fenwick tree of zeros.
fenny::update(&mut tree, 2 /* 0-based index */, 7);
fenny::update(&mut tree, 5, 20);
fenny::psum(&tree, 2);  // Returns 7.
fenny::psum(&tree, 5);  // Returns 27.

let root_p1 = fenny::get_root_p1(tree.len());
// Returns Some(5), as a[..=5].iter().sum() == 7 + 20 > 26.
fenny::first_larger(&tree, root_p1, 26);
fenny::first_larger(&tree, root_p1, 27); // Returns None.
```

## Slope-offset

What's better than a fenwick tree? Two of course! If we add another
tree into the mix, we can support an additional operation, the range
update. Also in the same juicy `O(log n)` time. The size of the range
doesn't matter, it can update the whole array in the same time as
updating a single element.

```rust
let mut slope = [0i64; 13];
let mut offset = [0i64; 13];
// Conceptually a[2..=5] += 3;
fenny::update_so(&mut slope, &mut offset, 2..=5, 3);
fenny::psum_so(&slope, &offset, 2); // Returns 3.
fenny::psum_so(&slope, &offset, 3); // Returns 6.
let root_p1 = fenny::get_root_p1(slope.len());
fenny::first_larger_so(&slope, &offset, root_p1, 3); // 3.
```

## Higher dimensions

This exists in 2d and 3d. `psum_2d(&tree, p)` computes a sum over all
points `q` satisfying `(q.y <= p.y && q.x <= p.x)`

```rust
use fenny::*;
let dim = Dim2{y: 5, x: 8};
let mut tree = vec![0; dim.y * dim.x];
update_2d(&mut tree, dim , Point2{y: 3, x: 6}, 7);
// Returns 0
psum_2d(&tree, dim, Point2{y: 2, x: 6});
// Returns 0
psum_2d(&tree, dim, Point2{y: 3, x: 5});
// Returns 7
psum_2d(&tree, dim, Point2{y: 3, x: 6});
// Returns 7
psum_2d(&tree, dim, Point2{y: 4, x: 7});
```
These operations run in `log(dim.x) * log (dim.y) * log (dim.z)`

## Higher dimension, slope offset, lexicographic

First a little background. This is something I came up with to help
with arithmetic coding. We have a probability density function over a
discretized space. Now to turn this into an efficient encoding, we
need to map our space to [0, 1), where every cell corresponds to a
segment with a length proportional to its probability. So lets imagine
we lay all the cells out in a line, and then we compute the prefix
sums of the probabilities. That will get us segments wohoo!

Ok, but one issue. We would like to efficiently manipulate the
probabilities of the cells, as they exist in 3d space. In particular
we would like to pick a box with corners (x<sub>0</sub>,
y<sub>0</sub>, z<sub>0</sub>) and (x<sub>1</sub>, y<sub>1</sub>,
z<sub>1</sub>), and increase the weight of all cells in that box in
one swoop. We could do range updates in 1d with a slope-offset
construction, anything similar for 3d? Yes there is! Turns out that by
having one slope per dimension, we can achieve our goal! This is
similar to an algorithm by [Mishra](https://arxiv.org/abs/1311.6093),
but we can get away a little cheaper by only suppporting lexicographic
queries rather than arbitrary boxes. During the decoding step it's
necessary to go from a number in [0, 1) and find the corresponding
cell. This can be done using `first_larger` function family.

With that background out of the way, let's look at the operations.

```rust
use fenny::*;
let dim = Dim3{z: 12, y: 6, x: 8};
let mut slope_z = vec![0i64; dim.z];
let mut slope_y = vec![0i64; dim.z * dim.y];
let mut slope_x = vec![0i64; dim.z * dim.y * dim.x];
let mut offset  = vec![0i64; dim.z * dim.y * dim.x];
let p0 = Point3{z: 3, y: 2, x: 1};
let p1 = Point3{z: 7, y: 5, x: 5};
let p2 = Point3{z: 4, y: 0, x: 0};
// Updates all points q with p0.x<=q.x<=p1.x && p0.y<=q.x<=p1.y && p0.z<=q.x<=p1.z
update_so_3d_lex(&mut slope_z, &mut slope_y, &mut slope_x, &mut offset, dim, p0, p1, 11);
// Our points out are laid out in lexicographic order according to z, y, x.
// Since 4, 0, 0 comes after (3, 2, 1), (3, 2, 2), ... (3, 5, 5) our result is 5 * 4 * 11;
psum_so_3d_lex(&slope_z, &slope_y, &slope_x, &offset, dim, p2);
// (4, 0, 0) has a psum of 5*4*11, but so does (3, 5, 5), so that's the result.
first_larger_so_3d_lex(&slope_z, &slope_y, &slope_x, &offset, dim, 5 * 4 * 11 - 1);
```

These operations run in `log(dim.x) * log (dim.y) * log (dim.z)`

## Acknowledgements

I found a topcoder
[article](https://www.topcoder.com/thrive/articles/Binary%20Indexed%20Trees)
written by boba5551 very helpful in understanding the Fenwick
trees. The diagrams were great.

The wikipedia
[article](https://en.wikipedia.org/w/index.php?title=Fenwick_tree) had
a good overview as well as useful information on how to use to do
range updates using two combined trees.

The [fenwick](https://crates.io/crates/fenwick) rust library is well
written and was an inspiration.
