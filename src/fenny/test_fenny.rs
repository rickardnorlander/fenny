mod fenny {
    #[cfg(test)]
    mod tests {
        use crate::fenny::*;
        use rand::Rng;
        use rand::SeedableRng;

        #[test]
        fn index_functions() {
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

            assert_eq!(get_root_p1(1), 1);
            assert_eq!(get_root_p1(2), 2);
            assert_eq!(get_root_p1(3), 2);
            assert_eq!(get_root_p1(4), 4);
            assert_eq!(get_root_p1(5), 4);
            assert_eq!(get_root_p1(6), 4);
            assert_eq!(get_root_p1(7), 4);
            assert_eq!(get_root_p1(8), 8);
        }

        #[test]
        #[should_panic]
        fn set_least_significant_zero_overflow() {
            set_least_significant_zero(usize::MAX);
        }

        #[test]
        #[should_panic]
        fn unset_trailing_ones_overflow() {
            unset_trailing_ones(usize::MAX);
        }

        #[test]
        #[should_panic]
        fn get_root_p1_underflow() {
            use crate::fenny::get_root_p1;
            get_root_p1(usize::MAX);
        }

        #[test]
        fn single_fenny() {
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
        fn first_larger_size_edgecases() {
            for size in [1, 2, 14, 15, 16, 16, 17] {
                let mut fenny_arr = vec![0; size];
                let rp1 = get_root_p1(size);
                fenny_arr[size - 1] = 1;
                assert_eq!(first_larger(&fenny_arr, rp1, -1), Some(0));
                assert_eq!(first_larger(&fenny_arr, rp1, 0), Some(size - 1));
                assert_eq!(first_larger(&fenny_arr, rp1, 1), None);
            }
        }

        #[test]
        fn single_fenny_random() {
            use rand::Rng;
            use rand::SeedableRng;
            let mut r = rand::rngs::StdRng::from_seed([0;32]);

            for _ in 0..100 {
                let s = r.gen_range(1..20);
                let rp1 = get_root_p1(s);
                let mut fenny_arr = vec![0; s];
                let mut plain_arr = vec![0; s];
                for _ in 0..10 {
                    let ind = r.gen_range(0..s);
                    let val = r.gen_range(0..=10);
                    update(&mut fenny_arr, ind, val);
                    plain_arr[ind] += val;
                }
                let mut acc = 0;
                for i in 0..s {
                    acc += plain_arr[i];
                    assert_eq!(psum(&fenny_arr, i), acc);
                    if plain_arr[i] != 0 {
                        assert_eq!(first_larger(&fenny_arr, rp1, acc-1), Some(i));
                    }
                }
            }
        }

        #[test]
        fn fenny_so() {
            let mut fenny_o = vec![0; 13];
            let mut fenny_s = vec![0; 13];

            let rp1 = get_root_p1(fenny_s.len());
            assert_eq!(rp1, 8);

            update_so(&mut fenny_s, &mut fenny_o, 3..=5, 7);
            assert_eq!(psum_so(&fenny_s, &fenny_o, 1), 0);
            assert_eq!(psum_so(&fenny_s, &fenny_o, 3), 7);
            assert_eq!(psum_so(&fenny_s, &fenny_o, 4), 14);
            assert_eq!(psum_so(&fenny_s, &fenny_o, 5), 21);
            assert_eq!(psum_so(&fenny_s, &fenny_o, 6), 21);


            assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, -1), Some(0));
            assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, 0), Some(3));
            assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, 11), Some(4));
            assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, 17), Some(5));
            assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, 21), None);
        }

        #[test]
        fn first_larger_size_edgecases_so() {
            for size in [1, 2, 14, 15, 16, 16, 17] {
                let mut fenny_s = vec![0i64; size];
                let mut fenny_o = vec![0; size];
                let rp1 = get_root_p1(size);
                update_so(&mut fenny_s, &mut fenny_o, 0..=size-1, 1);
                assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, -1), Some(0));
                assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, size as i64 -1), Some(size - 1));
                assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, size as i64), None);
            }
        }

        #[test]
        fn fenny_so_random() {
            let mut r = rand::rngs::StdRng::from_seed([0;32]);

            for _ in 0..100 {
                let s = r.gen_range(1..20);
                let rp1 = get_root_p1(s);
                let mut fenny_s = vec![0; s];
                let mut fenny_o = vec![0; s];
                let mut plain_arr = vec![0; s];
                for _ in 0..10 {
                    let ind1 = r.gen_range(0..s);
                    let ind2 = r.gen_range(ind1..s);
                    let val = r.gen_range(0..=10);
                    for j in ind1..=ind2 {
                        plain_arr[j] += val;
                    }
                    update_so(&mut fenny_s, &mut fenny_o, ind1..=ind2, val);
                }
                let mut acc = 0;
                for i in 0..s {
                    acc += plain_arr[i];
                    assert_eq!(psum_so(&fenny_s, &fenny_o, i), acc);
                    if plain_arr[i] != 0 {
                        assert_eq!(first_larger_so(&fenny_s, &fenny_o, rp1, acc-1), Some(i));
                    }
                }
            }
        }

        #[test]
        fn dim2_random() {
            let mut r = rand::rngs::StdRng::from_seed([0;32]);
            for _ in 0..100 {
                let dim = Dim2{x:r.gen_range(1..10), y:r.gen_range(1..10)};
                let mut fenny = vec![0; dim.x * dim.y];
                let mut arr = vec![0; dim.x * dim.y];
                for _ in 0..10 {
                    let p = Point2{x:r.gen_range(0..dim.x), y:r.gen_range(0..dim.y)};
                    let val = r.gen_range(0..=10);
                    update_2d(&mut fenny, dim, p, val);
                    arr[dim.x * p.y + p.x] += val
                }
                for y in 0..dim.y {
                    for x in 1..dim.x {
                        arr[y * dim.x + x] += arr[y * dim.x + x-1];
                    }
                }
                for y in 1..dim.y {
                    for x in 0..dim.x {
                        arr[y * dim.x + x] += arr[(y-1) * dim.x + x];
                        assert_eq!(psum_2d(&fenny, dim, Point2{x, y}), arr[y * dim.x + x]);
                    }
                }
            }
        }

        #[test]
        fn dim2_so() {
            let dim = Dim2{x: 8, y: 7};
            let size = dim.x * dim.y;
            let mut fenny_y = vec![0; dim.y];
            let mut fenny_x = vec![0; size];
            let mut fenny_o = vec![0; size];
            let p0 = Point2{y:2, x:2};
            let p1 = Point2{y:3, x:4};

            update_so_2d_lex(&mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, 7);

            let mut myarr = vec![0; size];
            for y in p0.y..=p1.y {
                for x in p0.x..=p1.x {
                    myarr[y * dim.x + x] = 7;
                }
            }
            let mut cumsum = 0;
            for y in 0..dim.y {
                for x in 0..dim.x {
                    let p = Point2{y, x};
                    cumsum += myarr[y * dim.x + x];
                    assert_eq!(psum_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, p), cumsum);
                    if myarr[y * dim.x + x] > 0 {
                        assert_eq!(first_larger_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, cumsum - 1), Some(p));
                    }
                }
            }
            assert_eq!(first_larger_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, cumsum), None);
        }
        #[test]
        fn dim2_so_random() {
            let mut r = rand::rngs::StdRng::from_seed([0;32]);
            for _ in 0..100 {
                let dim = Dim2{x: r.gen_range(1..20), y: r.gen_range(1..20)};
                let size = dim.x * dim.y;
                let mut fenny_y = vec![0; dim.y];
                let mut fenny_x = vec![0; size];
                let mut fenny_o = vec![0; size];
                let mut myarr = vec![0; size];

                for _ in 0..2 {
                    let p0 = Point2{y:r.gen_range(0..dim.y), x:r.gen_range(0..dim.x)};
                    let p1 = Point2{y:r.gen_range(p0.y..dim.y), x:r.gen_range(p0.x..dim.x)};
                    let val = r.gen_range(0..=10);

                    update_so_2d_lex(&mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, val);

                    for y in p0.y..=p1.y {
                        for x in p0.x..=p1.x {
                            myarr[y * dim.x + x] += val;
                        }
                    }
                }
                let mut cumsum = 0;
                for y in 0..dim.y {
                    for x in 0..dim.x {
                        cumsum += myarr[y * dim.x + x];
                        let p = Point2{y, x};
                        assert_eq!(psum_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, p), cumsum);
                        if myarr[y * dim.x + x] > 0 {
                            assert_eq!(first_larger_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, cumsum - 1), Some(p));
                        }
                    }
                }
                assert_eq!(first_larger_so_2d_lex(&fenny_y, &fenny_x, &fenny_o, dim, cumsum), None);
            }
        }

        #[test]
        fn dim3_so() {
            let test_cases = [[1, 1, 1, 0, 0, 0, 0, 0, 0],
                              [4, 1, 1, 0, 0, 0, 3, 0, 0],
                              [1, 4, 1, 0, 0, 0, 0, 3, 0],
                              [1, 1, 4, 0, 0, 0, 0, 0, 3],
                              [4, 4, 4, 0, 0, 0, 3, 3, 0],
                              [4, 4, 4, 0, 0, 0, 0, 3, 3],
                              [4, 4, 4, 0, 0, 0, 3, 0, 3],
                              [4, 4, 4, 1, 0, 0, 3, 3, 3],
                              [4, 4, 4, 0, 1, 0, 3, 3, 3],
                              [4, 4, 4, 0, 0, 1, 3, 3, 3],
                              [4, 4, 4, 1, 1, 0, 3, 3, 3],
                              [4, 4, 4, 0, 1, 1, 3, 3, 3],
                              [4, 4, 4, 1, 0, 1, 3, 3, 3],
            ];
            for tc in test_cases {
                let dim = Dim3{z: tc[0], y: tc[1], x: tc[2]};
                let p0 = Point3{z: tc[3], y: tc[4], x: tc[5]};
                let p1 = Point3{z: tc[6], y: tc[7], x: tc[8]};
                let size = dim.x * dim.y * dim.z;
                let mut fenny_z = vec![0; dim.z];
                let mut fenny_y = vec![0; dim.z * dim.y];
                let mut fenny_x = vec![0; dim.x * dim.y * dim.z];
                let mut fenny_o = vec![0; dim.x * dim.y * dim.z];
                let mut myarr = vec![0; size];

                update_so_3d_lex(&mut fenny_z, &mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, 3) ;
                for z in p0.z..=p1.z {
                    for y in p0.y..=p1.y {
                        for x in p0.x..=p1.x {
                            myarr[z * dim.y * dim.x + y * dim.x + x] += 3;
                        }
                    }
                }
                let mut cumsum = 0;
                for z in 0..dim.z {
                    for y in 0..dim.y {
                        for x in 0..dim.x {
                            let p = Point3{z, y, x};
                            let v = myarr[z * dim.y * dim.x + y * dim.x + x];
                            cumsum += v;
                            assert_eq!(psum_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, p), cumsum);
                            if v > 0 {
                                assert_eq!(first_larger_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, cumsum - 1), Some(p));
                            }
                        }
                    }
                }
                assert_eq!(first_larger_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, cumsum), None);
            }
        }

        #[test]
        fn dim3_so_random() {
            let mut r = rand::rngs::StdRng::from_seed([0;32]);
            for _ in 0..10 {
                let dim =  Dim3{z: r.gen_range(1..50), y: r.gen_range(1..50), x: r.gen_range(1..50)};
                let mut fenny_z = vec![0; dim.z];
                let mut fenny_y = vec![0; dim.z * dim.y];
                let mut fenny_x = vec![0; dim.x * dim.y * dim.z];
                let mut fenny_o = vec![0; dim.x * dim.y * dim.z];
                let mut myarr = vec![0; dim.x * dim.y * dim.z];

                for _ in 0..10 {
                    let p0 = Point3{z: r.gen_range(0..dim.z),
                                    y: r.gen_range(0..dim.y),
                                    x: r.gen_range(0..dim.x)};

                    let p1 = Point3{z: r.gen_range(p0.z..dim.z),
                                    y: r.gen_range(p0.y..dim.y),
                                    x: r.gen_range(p0.x..dim.x)};

                    let v = r.gen_range(0..10);
                    update_so_3d_lex(&mut fenny_z, &mut fenny_y, &mut fenny_x, &mut fenny_o, dim, p0, p1, v) ;
                    for z in p0.z..=p1.z {
                        for y in p0.y..=p1.y {
                            for x in p0.x..=p1.x {
                                myarr[z * dim.y * dim.x + y * dim.x + x] += v;
                            }
                        }
                    }
                }
                let mut cumsum = 0;
                for z in 0..dim.z {
                    for y in 0..dim.y {
                        for x in 0..dim.x {
                            let v = myarr[z * dim.y * dim.x + y * dim.x + x];
                            let p = Point3{z, y, x};
                            cumsum += v;
                            assert_eq!(psum_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, p), cumsum);
                            if v > 0 {
                                assert_eq!(first_larger_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, cumsum - 1), Some(p));
                            }
                        }
                    }
                }
                assert_eq!(first_larger_so_3d_lex(&fenny_z, &fenny_y, &fenny_x, &fenny_o, dim, cumsum), None);
            }
        }
    }
}
