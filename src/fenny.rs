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
    let shift = 63 - (size+1).leading_zeros();
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

#[derive(Debug, Copy, Clone)]
pub struct Point2 {
    pub x: usize,
    pub y: usize,
}

pub fn psum_2d<T: Value>(fenny: &[T], dim: Dim2, p: Point2) -> T {
    let mut ret = T::default();
    for y2 in query_inds(p.y) {
        for x2 in query_inds(p.x) {
            ret += fenny[y2 * dim.x + x2];
        }
    }
    return ret;
}

pub fn update_2d<T: Value>(fenny: &mut [T], dim: Dim2,  p: Point2, val: T) {
    for y2 in update_inds(p.y, dim.y) {
        for x2 in update_inds(p.x, dim.x) {
            fenny[y2 * dim.x + x2] += val;
        }
    }
}



// 2d slope offset trees.

pub fn so_psum_2d_linear<T: Value>(f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim2, p: Point2) -> T {
    let s_y = psum_2d(f_slope_y, dim, p);
    let s_x = psum_2d(f_slope_x, dim, p);
    let o = psum_2d(f_offset, dim, p);

    return s_y * T::from_usize(p.y) + s_x * T::from_usize(p.x) + o;
}

fn helper_2d<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim2, p: Point2, val: T, ind_: usize, inclusive: bool)  {
    let ind = T::from_usize(ind_);
    let offset = if inclusive {
        val - val * ind
    } else {
        -val * ind
    };
    update_2d(f_slope, dim, p, val);
    update_2d(f_offset, dim, p, offset);
}

pub fn so_update_2d_linear<T: Value>(f_slope_y: &mut [T], f_slope_x: &mut [T], f_offset: &mut [T], dim: Dim2,
                           p0: Point2, p1: Point2, val: T) {
    assert!(p0.y <= p1.y && p1.y < dim.y);
    assert!(p0.x <= p1.x && p1.x < dim.x);
    let strip_value = val * T::from_usize(p1.x - p0.x + 1);

    helper_2d(f_slope_y, f_offset, dim, Point2{y:p0.y, x:0}, strip_value, p0.y, false);
    if p1.y + 1 < dim.y {
        helper_2d(f_slope_y, f_offset, dim, Point2{y:p1.y+1, x:0}, -strip_value, p1.y+1, false);
    }

    helper_2d(f_slope_x, f_offset, dim, Point2{y:p0.y, x:p0.x}, val, p0.x, true);
    helper_2d(f_slope_x, f_offset, dim, Point2{y:p0.y, x:p1.x}, -val, p1.x, false);
    if p1.y + 1 < dim.y {
        helper_2d(f_slope_x, f_offset, dim, Point2{y:p1.y+1, x:p0.x}, -val, p0.x, true);
        helper_2d(f_slope_x, f_offset, dim, Point2{y:p1.y+1, x:p1.x}, val, p1.x, false);
    }
}



// 3d Fenwick trees.

#[derive(Debug, Copy, Clone)]
pub struct Dim3 {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}

#[derive(Debug, Copy, Clone)]
pub struct Point3 {
    pub x: usize,
    pub y: usize,
    pub z: usize,
}


pub fn update_3d<T: Value>(fenny: &mut [T], dim: Dim3, p: Point3, val: T) {
    for z in update_inds(p.z, dim.z) {
        for y in update_inds(p.y, dim.y) {
            for x in update_inds(p.x, dim.x) {
                fenny[z * dim.y * dim.x + y * dim.x + x] += val;
            }
        }
    }
}

pub fn psum_3d<T: Value>(fenny: &[T], dim: Dim3, p: Point3) -> T {
    let mut ret = T::default();
    for z in query_inds(p.z) {
        for y in query_inds(p.y) {
            for x in query_inds(p.x) {
                ret += fenny[z * dim.y * dim.x + y * dim.x + x];
            }
        }
    }
    return ret;
}



// 3d slope offset trees.

pub fn so_psum_3d_linear<T: Value>(f_slope_z: &[T], f_slope_y: &[T], f_slope_x: &[T], f_offset: &[T], dim: Dim3, p: Point3) -> T {
    let s_z = psum_3d(f_slope_z, dim, p);
    let s_y = psum_3d(f_slope_y, dim, p);
    let s_x = psum_3d(f_slope_x, dim, p);
    let o = psum_3d(f_offset, dim, p);

    return s_z * T::from_usize(p.z) + s_y * T::from_usize(p.y) + s_x * T::from_usize(p.x) + o;
}

pub fn helper<T: Value>(f_slope: &mut [T], f_offset: &mut [T], dim: Dim3, p: Point3, val: T, ind_: usize, inclusive: bool) {
    let ind = T::from_usize(ind_);
    let offset = if inclusive {
        val - val * ind
    } else {
        -val * ind
    };
    update_3d(f_slope, dim, p, val);
    update_3d(f_offset, dim, p, offset);
}

pub fn so_update_3d_linear<T: Value>(f_slope_z: &mut [T], f_slope_y: &mut [T], f_slope_x: &mut [T], f_offset: &mut [T], dim: Dim3, p0: Point3, p1: Point3, val: T) {
    assert!(p0.z <= p1.z && p1.z < dim.z);
    assert!(p0.y <= p1.y && p1.y < dim.y);
    assert!(p0.x <= p1.x && p1.x < dim.x);

    // Contribution from complete xy-slabs.
    let slab_value = T::from_usize(p1.x - p0.x + 1) * T::from_usize(p1.y - p0.y + 1) * val;
    helper(f_slope_z, f_offset, dim, Point3{z:p0.z, y:0, x:0}, slab_value, p0.z, false);
    if p1.z + 1 < dim.z {
        helper(f_slope_z, f_offset, dim, Point3{z:p1.z+1, y:0, x:0}, -slab_value, p1.z+1, false);
    }

    // Contribution from complete x-strips.
    let strip_value = T::from_usize(p1.x - p0.x + 1) * val;
    helper(f_slope_y, f_offset, dim, Point3{z:p0.z, y:p0.y, x:0}, strip_value, p0.y, false);
    if p1.y + 1 < dim.y {
        helper(f_slope_y, f_offset, dim, Point3{z:p0.z, y:p1.y+1, x:0}, -strip_value, p1.y+1, false);
    }
    if p1.z + 1 < dim.z {
        helper(f_slope_y, f_offset, dim, Point3{z:p1.z+1, y:p0.y, x:0}, -strip_value, p0.y, false);
        if p1.y + 1 < dim.y {
            helper(f_slope_y, f_offset, dim, Point3{z:p1.z+1, y:p1.y+1, x:0}, strip_value, p1.y+1, false);
        }
    }

    // Contribution from x
    helper(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p0.y, x:p0.x},  val, p0.x, true);
    helper(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p0.y, x:p1.x}, -val, p1.x, false);
    if p1.y + 1 < dim.y {
        helper(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p1.y + 1, x:p0.x},  -val, p0.x, true);
        helper(f_slope_x, f_offset, dim, Point3{z: p0.z, y:p1.y + 1, x:p1.x},   val, p1.x, false);
    }
    if p1.z + 1 < dim.z {
        helper(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p0.y, x:p0.x},  -val, p0.x, true);
        helper(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p0.y, x:p1.x},   val, p1.x, false);
        if p1.y + 1 < dim.y {
            helper(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p1.y + 1, x:p0.x},   val, p0.x, true);
            helper(f_slope_x, f_offset, dim, Point3{z: p1.z+1, y:p1.y + 1, x:p1.x},  -val, p1.x, false);
        }
    }
}



// Tests

#[cfg(test)]
mod tests {
    #[test]
    fn index_functions() {
        use crate::fenny::set_least_significant_zero;
        use crate::fenny::unset_trailing_ones;

        assert_eq!(set_least_significant_zero(0), 1);
        assert_eq!(set_least_significant_zero(0b1100), 0b1101);
        assert_eq!(set_least_significant_zero(0b1101), 0b1111);
        assert_eq!(set_least_significant_zero(0b1111), 0b11111);

        assert_eq!(unset_trailing_ones(1), 0);
        assert_eq!(unset_trailing_ones(0), 0);
        assert_eq!(unset_trailing_ones(0b1000), 0b1000);
        assert_eq!(unset_trailing_ones(0b1001), 0b1000);
        assert_eq!(unset_trailing_ones(0b1011), 0b1000);
        assert_eq!(unset_trailing_ones(0b0011), 0);
    }
    #[test]
    fn single_fenny() {
        use crate::fenny::update;
        use crate::fenny::psum;
        use crate::fenny::first_larger;
        use crate::fenny::get_root_p1;

        let mut fenny_arr = vec![0; 13];

        let rp1 = get_root_p1(fenny_arr.len());
        assert_eq!(rp1, 8);

        update(&mut fenny_arr, 3, 17);
        assert_eq!(psum(&mut fenny_arr, 2), 0);
        assert_eq!(psum(&mut fenny_arr, 3), 17);
        assert_eq!(psum(&mut fenny_arr, 4), 17);

        assert_eq!(first_larger(&fenny_arr, rp1, -1), Some(0));
        assert_eq!(first_larger(&fenny_arr, rp1, 0), Some(3));
        assert_eq!(first_larger(&fenny_arr, rp1, 16), Some(3));
        assert_eq!(first_larger(&fenny_arr, rp1, 17), None);

        update(&mut fenny_arr, 4, 2);
        assert_eq!(psum(&mut fenny_arr, 3), 17);
        assert_eq!(psum(&mut fenny_arr, 4), 19);
        assert_eq!(psum(&mut fenny_arr, 5), 19);
        assert_eq!(first_larger(&fenny_arr, rp1, 18), Some(4));
        assert_eq!(first_larger(&fenny_arr, rp1, 19), None);
    }
    #[test]
    fn fenny_so() {
        use crate::fenny::so_update;
        use crate::fenny::so_psum;
        use crate::fenny::so_first_larger;
        use crate::fenny::get_root_p1;

        let mut fenny_o = vec![0; 13];
        let mut fenny_s = vec![0; 13];

        let rp1 = get_root_p1(fenny_s.len());
        assert_eq!(rp1, 8);

        so_update(&mut fenny_s, &mut fenny_o, 3..=5, 7);
        assert_eq!(so_psum(&fenny_s, &fenny_o, 1), 0);
        assert_eq!(so_psum(&fenny_s, &fenny_o, 3), 7);
        assert_eq!(so_psum(&fenny_s, &fenny_o, 4), 14);
        assert_eq!(so_psum(&fenny_s, &fenny_o, 5), 21);
        assert_eq!(so_psum(&fenny_s, &fenny_o, 6), 21);


        assert_eq!(so_first_larger(&fenny_s, &fenny_o, rp1, -1), Some(0));
        assert_eq!(so_first_larger(&fenny_s, &fenny_o, rp1, 0), Some(3));
        assert_eq!(so_first_larger(&fenny_s, &fenny_o, rp1, 11), Some(4));
        assert_eq!(so_first_larger(&fenny_s, &fenny_o, rp1, 17), Some(5));
        assert_eq!(so_first_larger(&fenny_s, &fenny_o, rp1, 21), None);
    }
    #[test]
    fn dim2() {
        use crate::fenny::so_update_2d_linear;
        use crate::fenny::so_psum_2d_linear;
        use crate::fenny::Point2;
        use crate::fenny::Dim2;

        let dim = Dim2{x: 8, y: 7};
        let size = dim.x * dim.y;
        let mut fenny_y = vec![0; size];
        let mut fenny_x = vec![0; size];
        let mut fenny_o = vec![0; size];
        let p0 = Point2{y:2, x:2};
        let p1 = Point2{y:3, x:4};

        so_update_2d_linear(&mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, 7);

        let mut myarr = vec![0; size];
        for y in p0.y..=p1.y {
            for x in p0.x..=p1.x {
                myarr[y * dim.x + x] = 7;
            }
        }
        let mut cumsum = 0;
        for y in 0..dim.y {
            for x in 0..dim.x {
                cumsum += myarr[y * dim.x + x];
                //println!("{} {} {}", y, x, cumsum);
                assert_eq!(so_psum_2d_linear(&fenny_y, &fenny_x, &fenny_o, dim, Point2{y, x}), cumsum);
            }
        }
    }
    #[test]
    fn dim3() {
        use crate::fenny::so_update_3d_linear;
        use crate::fenny::so_psum_3d_linear;
        use crate::fenny::{Dim3,Point3};

        let dim =  Dim3{z:11, y: 11, x: 7};
        let size = dim.x * dim.y * dim.z;
        let mut fenny_z = vec![0; size];
        let mut fenny_y = vec![0; size];
        let mut fenny_x = vec![0; size];
        let mut fenny_o = vec![0; size];

        let p0 = Point3{z:4,y:3,x:2};
        let p1 = Point3{z:6,y:8,x:5};

        so_update_3d_linear(&mut fenny_z, &mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, 7);

        let mut myarr = vec![0; size];
        for z in p0.z..=p1.z {
            for y in p0.y..=p1.y {
                for x in p0.x..=p1.x {
                    myarr[z * dim.y * dim.x + y * dim.x + x] = 7;
                }
            }
        }
        let mut cumsum = 0;
        for z in 0..dim.z {
            for y in 0..dim.y {
                for x in 0..dim.x {
                    cumsum += myarr[z * dim.y * dim.x + y * dim.x + x];
                    let p = Point3{z, y, x};
                    println!("{} {} {} {}", z, y, x, cumsum);
                    assert_eq!(so_psum_3d_linear(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, p), cumsum);
                }
            }
        }
    }
}

use std::num::Wrapping;

macro_rules! impl_primitive {
    ($primitive:ty) => {
        impl Value for $primitive {
            fn from_usize(v: usize) -> Self {
                v as Self
            }
        }
    };
}
macro_rules! impl_wrapping {
    ($primitive:ty) => {
        impl Value for Wrapping<$primitive> {
            fn from_usize(v: usize) -> Self {
                Wrapping(v as $primitive)
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
impl_wrapping!(i16);
impl_wrapping!(i32);
impl_wrapping!(i64);
impl_wrapping!(i128);
impl_wrapping!(u16);
impl_wrapping!(u32);
impl_wrapping!(u64);
impl_wrapping!(u128);
