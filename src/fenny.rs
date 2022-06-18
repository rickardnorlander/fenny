// BSD 2-Clause License
//
// Copyright (c) 2022, Rickard Norlander
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//
//
//  Conventions:
//
//  There are a few variants of Fenwick trees. This uses 0-based
//  indexing, meaning that
//                                  tree[10100111]
//
// stores sum of all numbers of the form 10100xxx
// Basically the trailing one-bits act like wildcards.
//
// Examples
// So if we want to modify the number with index
//               10100101
// We have to update any tree node that sums over that, namely
//          tree[10100101]
//          tree[10100111]
//          tree[10101111]
//          tree[10111111]
//          tree[11111111]
//
// If we want to query prefix sum of
//               01011010
// Then note that
//           pre(01011010)=
//           sum(01011010)+pre(01011001)
//           sum(01011010)+sum(0101100x)+pre(01010111)
//           sum(01011010)+sum(0101100x)+sum(01010xxx)+pre(01001111)
//           sum(01011010)+sum(0101100x)+sum(01010xxx)+sum(0100xxxx)+pre(00111111)
//           sum(01011010)+sum(0101100x)+sum(01010xxx)+sum(0100xxxx)+sum(00xxxxxx)
//
// Each of which correspond to some node in the tree.


use std::num::Wrapping;



// Value type

pub trait Value :
    Copy
    + std::default::Default
    + std::ops::AddAssign
    + std::ops::Add<Output=Self>
    + std::iter::Sum
    + std::ops::Neg<Output=Self>
    + std::ops::Sub<Output=Self>
    + std::ops::Mul<Output=Self>
    + std::cmp::PartialOrd
{
    fn from_usize(v: usize) -> Self;
}
macro_rules! impl_primitive {
    ($primitive:ty) => {
        impl Value for $primitive {
            fn from_usize(v: usize) -> Self {
                v as Self
            }
        }
    };
}
impl_primitive!(i16);
impl_primitive!(i32);
impl_primitive!(i64);
impl_primitive!(i128);
impl_primitive!(f32);
impl_primitive!(f64);
macro_rules! impl_wrapping {
    ($primitive:ty) => {
        impl Value for Wrapping<$primitive> {
            fn from_usize(v: usize) -> Self {
                Wrapping(v as $primitive)
            }
        }
    };
}
impl_wrapping!(i16);
impl_wrapping!(i32);
impl_wrapping!(i64);
impl_wrapping!(i128);
impl_wrapping!(u16);
impl_wrapping!(u32);
impl_wrapping!(u64);
impl_wrapping!(u128);



// Indexing utilities

fn set_least_significant_zero(ind: usize) -> usize {
    return ind | (ind + 1);
}

fn unset_trailing_ones(ind: usize) -> usize {
    return ind & (ind + 1);
}

pub const fn get_root_p1(size: usize) -> usize {
    if size == 0 {
        panic!();
    }
    let shift = 63 - size.leading_zeros();
    (1 << shift) as usize
}

fn query_inds(ind: usize) -> impl Iterator<Item=usize> {
    return std::iter::successors(Some(ind), |n| {
        let m = unset_trailing_ones(*n);
        if m > 0 {Some(m-1)} else {None}
    });
}



// Basic Fenwick tree.

pub fn psum<T: Value>(fenny: &[T], ind: usize) -> T {
    return query_inds(ind).map(|ind| fenny[ind]).sum();
}

fn update_inds(ind: usize, fenny_len: usize) -> impl Iterator<Item=usize>  {
    return std::iter::successors(Some(ind), move |n| {
        let m = set_least_significant_zero(*n);
        if m < fenny_len {Some(m)} else {None}
    });
}

pub fn update<T: Value>(fenny: &mut [T], ind: usize, val: T) {
    for ind in update_inds(ind, fenny.len()) {
        fenny[ind] += val;
    }
}

struct BSInds {
    size: usize,
    low_p1: usize,
    step: usize,
}

impl BSInds {
    pub fn from_root (root_p1: usize, size: usize) -> Self {
        Self {size, step: root_p1, low_p1: 0}
    }
    pub fn next(&mut self) -> Option<usize> {
        while self.step > 0 {
            let next = self.low_p1 + self.step - 1;
            self.step /= 2;
            if next < self.size {
                return Some(next);
            }
        }
        return None;
    }
    pub fn higher(&mut self) {
        self.low_p1 += self.step * 2;
    }
}

// Finds the smallest ind so that val < arr[i]
pub fn first_larger<T: Value>(fenny: &[T], root_p1: usize, val: T) -> Option<usize> {
    let mut inds = BSInds::from_root(root_p1, fenny.len());
    let mut result = None;
    let mut low_val = T::default();
    while let Some(ind) = inds.next() {
        let possible = fenny[ind] + low_val;
        if possible <= val {
            low_val = possible;
            inds.higher();
        } else {
            result = Some(ind);
        }
    }
    result
}



// Slope-offset Fenwick tree.

pub fn so_update<T: Value>(slope: &mut [T], offset: &mut [T], range: std::ops::RangeInclusive<usize>, val: T) {
    let (start, end) = (*range.start(), *range.end());
    update(slope, start, val);
    update(slope, end, -val);
    update(offset, start,  val - val * T::from_usize(start));
    update(offset, end,  val * T::from_usize(end));
}

pub fn so_psum<T: Value>(slope: &[T], offset: &[T], ind: usize) -> T {
    return psum(slope, ind) * T::from_usize(ind) +
           psum(offset, ind);
}

// Goal, find the smallest ind so that val < arr[i]
pub fn so_first_larger<T: Value>(f_slope: &[T], f_offset: &[T], root_p1: usize, val: T) -> Option<usize> {
    let mut inds = BSInds::from_root(root_p1, f_slope.len());
    let mut low_offset = T::default();
    let mut low_slope = T::default();
    let mut result = None;
    while let Some(ind) = inds.next() {
        let possible_offset = f_offset[ind] + low_offset;
        let possible_slope =  f_slope [ind] + low_slope;
        let possible = possible_slope * T::from_usize(ind) + possible_offset;
        if possible <= val {
            low_offset = possible_offset;
            low_slope = possible_slope;
            inds.higher();
        } else {
            result = Some(ind);
        }
    }
    result
}



// 2d Fenwick trees.

#[derive(Debug, Copy, Clone)]
pub struct Dim2 {
    pub x: usize,
    pub y: usize,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point2 {
    pub x: usize,
    pub y: usize,
}

pub fn psum_2d<T: Value>(fenny: &[T], dim: Dim2, p: Point2) -> T {
    let mut ret = T::default();
    for x2 in query_inds(p.x) {
        for y2 in query_inds(p.y) {
            ret += fenny[x2 * dim.y + y2];
        }
    }
    return ret;
}

pub fn update_2d<T: Value>(fenny: &mut [T], dim: Dim2,  p: Point2, val: T) {
    for x2 in update_inds(p.x, dim.x) {
        for y2 in update_inds(p.y, dim.y) {
            fenny[x2 * dim.y + y2] += val;
        }
    }
}



// 2d slope offset trees.

pub fn so_psum_2d_lex<T: Value>(f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim2, p: Point2) -> T {
    let s_y = psum(f_slope_y, p.y);
    let s_x = psum_2d(f_slope_x, dim, p);
    let o = psum_2d(f_offset, dim, p);

    return s_y * T::from_usize(p.y) + s_x * T::from_usize(p.x) + o;
}

fn helper_2y<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim2, y: usize, val: T)  {
    let ind = T::from_usize(y);
    let offset = -val * ind;
    update(f_slope, y, val);
    update_2d(f_offset, dim, Point2{y, x:0}, offset);
}
fn helper_2x<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim2, p: Point2, val: T, ind_: usize, inclusive: bool)  {
    let ind = T::from_usize(ind_);
    let offset = if inclusive {
        val - val * ind
    } else {
        -val * ind
    };
    update_2d(f_slope, dim, p, val);
    update_2d(f_offset, dim, p, offset);
}

pub fn so_update_2d_lex<T: Value>(f_slope_y: &mut [T], f_slope_x: &mut [T], f_offset: &mut [T], dim: Dim2,
                           p0: Point2, p1: Point2, val: T) {
    assert!(p0.y <= p1.y && p1.y < dim.y);
    assert!(p0.x <= p1.x && p1.x < dim.x);
    let strip_value = val * T::from_usize(p1.x - p0.x + 1);

    helper_2y(f_slope_y, f_offset, dim, p0.y, strip_value);
    if p1.y + 1 < dim.y {
        helper_2y(f_slope_y, f_offset, dim, p1.y+1, -strip_value);
    }

    helper_2x(f_slope_x, f_offset, dim, Point2{y:p0.y, x:p0.x}, val, p0.x, true);
    helper_2x(f_slope_x, f_offset, dim, Point2{y:p0.y, x:p1.x}, -val, p1.x, false);
    if p1.y + 1 < dim.y {
        helper_2x(f_slope_x, f_offset, dim, Point2{y:p1.y+1, x:p0.x}, -val, p0.x, true);
        helper_2x(f_slope_x, f_offset, dim, Point2{y:p1.y+1, x:p1.x}, val, p1.x, false);
    }
}

fn marginalize_out_y<T: Value>(fenny: &[T], dim: Dim2, y: usize, x: usize) -> T {
    return psum(&fenny[x*dim.y..], y);
}

pub fn so_first_larger_2d_lex<T: Value>(f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim2, val: T) -> Option<Point2> {
    let root_y_p1 = get_root_p1(dim.y);
    let mut yinds = BSInds::from_root(root_y_p1, dim.y);
    let mut result = None;

    let mut low_offset = T::default();
    let mut low_slope_y = T::default();

    while let Some(y) = yinds.next() {
        let possible_slope_y = low_slope_y + f_slope_y[y];
        let possible_offset = low_offset + f_offset[y];
        let possible = possible_slope_y * T::from_usize(y) + possible_offset;
        if possible <= val {
            low_slope_y = possible_slope_y;
            low_offset = possible_offset;
            yinds.higher();
        } else {
            result = Some(Point2{y, x:0});
        }
    }
    if result == Some(Point2{y: 0, x: 0}) {
        return result;
    }
    // We know now that prefix_sum[y-1] <= val < prefix_sum[y]
    // Our answer could either be (y, 0) or it could be (y-1, something)

    let root_x_p1 = get_root_p1(dim.x);
    let mut xinds = BSInds::from_root(root_x_p1, dim.x);

    let y = result.map(|yx| yx.y).unwrap_or(dim.y) - 1;
    let ypart = T::from_usize(y) * psum(f_slope_y, y);

    low_offset = T::default();
    let mut low_slope_x = T::default();

    while let Some(x) = xinds.next() {
        let possible_slope_x = low_slope_x + marginalize_out_y(f_slope_x, dim, y, x);
        let possible_offset = low_offset + marginalize_out_y(f_offset, dim, y, x);
        let possible = ypart + possible_slope_x * T::from_usize(x) + possible_offset;
        if possible <= val {
            low_offset = possible_offset;
            low_slope_x = possible_slope_x;
            xinds.higher();
        } else {
            result = Some(Point2{y, x})
        }
    }
    result
}


// 3d Fenwick trees.

#[derive(Debug, Copy, Clone)]
pub struct Dim3 {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point3 {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}


pub fn update_3d<T: Value>(fenny: &mut [T], dim: Dim3, p: Point3, val: T) {
    for x in update_inds(p.x, dim.x) {
        for y in update_inds(p.y, dim.y) {
            for z in update_inds(p.z, dim.z) {
                fenny[z + y * dim.z + x * dim.y * dim.z] += val;
            }
        }
    }
}

pub fn psum_3d<T: Value>(fenny: &[T], dim: Dim3, p: Point3) -> T {
    let mut ret = T::default();
    for x in query_inds(p.x) {
        for y in query_inds(p.y) {
            for z in query_inds(p.z) {
                ret += fenny[z + y * dim.z + x * dim.y * dim.z];
            }
        }
    }
    return ret;
}



// 3d slope offset trees.

pub fn so_psum_3d_lex<T: Value>(f_slope_z: &[T], f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim3, p: Point3) -> T {
    let s_z = psum(f_slope_z, p.z);
    let s_y = psum_2d(f_slope_y, Dim2{y:dim.z,x:dim.y}, Point2{y:p.z,x:p.y});
    let s_x = psum_3d(f_slope_x, dim, p);
    let o = psum_3d(f_offset, dim, p);

    return s_z * T::from_usize(p.z) + s_y * T::from_usize(p.y) + s_x * T::from_usize(p.x) + o;
}

fn helper_3z<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim3, z: usize, val: T) {
    let ind = T::from_usize(z);
    let offset = -val * ind;
    update(f_slope, z, val);
    update_3d(f_offset, dim, Point3{z:z, y:0, x:0}, offset);
}

fn helper_3y<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim3, p: Point3, val: T) {
    let ind = T::from_usize(p.y);
    let offset = -val * ind;
    update_2d(f_slope, Dim2{y:dim.z, x: dim.y}, Point2{y:p.z, x:p.y}, val);
    update_3d(f_offset, dim, Point3{z:p.z, y:p.y, x:0}, offset);
}

fn helper_3x<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim3, p: Point3, val: T, inclusive: bool) {
    let ind = T::from_usize(p.x);
    let offset = if inclusive {
        val - val * ind
    } else {
        -val * ind
    };
    update_3d(f_slope, dim, p, val);
    update_3d(f_offset, dim, p, offset);
}

pub fn so_update_3d_lex<T: Value>(f_slope_z: &mut [T], f_slope_y: &mut [T], f_slope_x: &mut [T], f_offset: &mut [T], dim: Dim3, p0: Point3, p1: Point3, val: T) {
    assert!(p0.z <= p1.z && p1.z < dim.z);
    assert!(p0.y <= p1.y && p1.y < dim.y);
    assert!(p0.x <= p1.x && p1.x < dim.x);

    // Contribution from complete xy-slabs.
    let slab_value = T::from_usize(p1.x - p0.x + 1) * T::from_usize(p1.y - p0.y + 1) * val;
    helper_3z(f_slope_z, f_offset, dim, p0.z, slab_value);
    if p1.z + 1 < dim.z {
        helper_3z(f_slope_z, f_offset, dim, p1.z+1, -slab_value);
    }

    // Contribution from complete x-strips.
    let strip_value = T::from_usize(p1.x - p0.x + 1) * val;
    helper_3y(f_slope_y, f_offset, dim, Point3{z:p0.z, y:p0.y, x:0}, strip_value);
    if p1.y + 1 < dim.y {
        helper_3y(f_slope_y, f_offset, dim, Point3{z:p0.z, y:p1.y+1, x:0}, -strip_value);
    }
    if p1.z + 1 < dim.z {
        helper_3y(f_slope_y, f_offset, dim, Point3{z:p1.z+1, y:p0.y, x:0}, -strip_value);
        if p1.y + 1 < dim.y {
            helper_3y(f_slope_y, f_offset, dim, Point3{z:p1.z+1, y:p1.y+1, x:0}, strip_value);
        }
    }

    // Contribution from x
    helper_3x(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p0.y, x:p0.x},  val, true);
    helper_3x(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p0.y, x:p1.x}, -val, false);
    if p1.y + 1 < dim.y {
        helper_3x(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p1.y + 1, x:p0.x},  -val, true);
        helper_3x(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p1.y + 1, x:p1.x},   val, false);
    }
    if p1.z + 1 < dim.z {
        helper_3x(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p0.y, x:p0.x},  -val, true);
        helper_3x(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p0.y, x:p1.x},   val, false);
        if p1.y + 1 < dim.y {
            helper_3x(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p1.y + 1, x:p0.x},   val, true);
            helper_3x(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p1.y + 1, x:p1.x},  -val, false);
        }
    }
}

fn marginalize_out_z<T: Value>(fenny: &[T], dim: Dim3, z: usize, y: usize) -> T {
    return psum(&fenny[y*dim.z..], z);
}

fn marginalize_out_zy<T: Value>(fenny: &[T], dim: Dim3, z: usize, y: usize, x: usize) -> T{
    return psum_2d(&fenny[x * dim.z * dim.y..], Dim2{y: dim.z, x: dim.y}, Point2{y: z, x: y});
}


pub fn so_first_larger_3d_lex<T: Value>(f_slope_z: &[T], f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim3, val: T) -> Option<Point3> {
    let root_z_p1 = get_root_p1(dim.z);
    let mut zinds = BSInds::from_root(root_z_p1, dim.z);
    let mut result = None;
    let mut low_offset = T::default();
    let mut low_slope_z = T::default();

    while let Some(z) = zinds.next() {
        let possible_slope_z = low_slope_z + f_slope_z[z];
        let possible_offset = low_offset + f_offset[z];
        let possible = possible_slope_z * T::from_usize(z) + possible_offset;
        if possible <= val {
            low_slope_z = possible_slope_z;
            low_offset = possible_offset;
            zinds.higher();
        } else {
            result = Some(Point3{z, y:0, x:0});
        }
    }
    if result == Some(Point3{z:0, y:0, x:0}) {
        return result;
    }
    // We know now that prefix_sum[(z-1) * dim.y * dim.x] <= val < prefix_sum[z * dim.y * dim.x]
    // Our answer could either be (z, 0, 0) or it could be (z-1, something, something)
    let root_y_p1 = get_root_p1(dim.y);
    let mut yinds = BSInds::from_root(root_y_p1, dim.y);

    let z = result.map(|p| p.z).unwrap_or(dim.z) - 1;
    let zpart = T::from_usize(z) * psum(f_slope_z, z);

    low_offset = T::default();
    let mut low_slope_y = T::default();

    while let Some(y) = yinds.next() {
        let possible_slope_y = low_slope_y + marginalize_out_z(f_slope_y, dim, z, y);
        let possible_offset = low_offset + marginalize_out_z(f_offset, dim, z, y);
        let possible = zpart + possible_slope_y * T::from_usize(y) + possible_offset;
        if possible <= val {
            low_offset = possible_offset;
            low_slope_y = possible_slope_y;
            yinds.higher();
        } else {
            result = Some(Point3{z, y, x:0})
        }
    }
    if result == Some(Point3{z, y: 0, x: 0}) {
        return result;
    }
    let y = if result == None {
        dim.y - 1
    } else if result == Some(Point3{z: z+1, y: 0, x: 0}) {
        dim.y - 1
    } else {
        result.unwrap().y - 1
    };

    let root_x_p1 = get_root_p1(dim.x);
    let mut xinds = BSInds::from_root(root_x_p1, dim.x);

    let ypart = T::from_usize(y) * psum_2d(f_slope_y, Dim2{y: dim.z, x: dim.y}, Point2{y: z, x: y});

    low_offset = T::default();
    let mut low_slope_x = T::default();

    while let Some(x) = xinds.next() {
        let possible_slope_x = low_slope_x + marginalize_out_zy(f_slope_x, dim, z, y, x);
        let possible_offset = low_offset + marginalize_out_zy(f_offset, dim, z, y, x);
        let possible = zpart + ypart + possible_slope_x * T::from_usize(x) + possible_offset;
        if possible <= val {
            low_offset = possible_offset;
            low_slope_x = possible_slope_x;
            xinds.higher();
        } else {
            result = Some(Point3{z, y, x})
        }
    }
    result
}


mod test_fenny;
