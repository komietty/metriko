//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Port of libc++ 20.1 <__algorithm/sort.h>, <partial_sort.h>, <make_heap.h>, <sort_heap.h>, <sift_down.h>,
//--- <push_heap.h>, <pop_heap.h>, <max_element.h>, <min_element.h>, <minmax.h> and <queue> (priority_queue)
//--- (Apache-2.0 WITH LLVM-exception), as shipped with Apple clang 17 (MacOSX SDK, _LIBCPP_VERSION 200100).
//
// metriko's results depend on how libc++ orders elements that compare equal (or on comparators that are not strict
// weak orders). These functions perform the same comparisons in the same order and the same element moves as
// libc++, so they reproduce its output exactly.
//
// Scope: the reference build uses the default hardening mode (_LIBCPP_HARDENING_MODE_DEFAULT = NONE) and does not
// define _LIBCPP_DEBUG_RANDOMIZE_UNSPECIFIED_STABILITY, so neither the debug comparator wrapper (__debug_less) nor
// the range randomization is part of the behavior. `sort_by` is the path libc++ takes for a custom comparator or a
// projection (std::sort with a lambda, std::ranges::sort with a lambda or a projection): __introsort with
// _UseBitSetPartition = false, i.e. __partition_with_equals_on_right. The branchless / bitset path is only taken for
// arithmetic element types with std::less / std::greater, where equal elements are indistinguishable.
//
// libc++ moves values through a temporary (the pivot, the inserted element, the heap top). Here the value stays in
// the slice and travels by swaps / rotations instead: every comparison sees the same value as in libc++, and the
// final arrangement is the same.

/// `std::sort(first, last, less)` / `std::ranges::sort(r, less)` / `std::ranges::sort(r, {}, proj)` (with
/// `less(a, b) = proj(a) < proj(b)`).
pub fn sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    let n = v.len();
    // __sort_dispatch: __depth_limit = 2 * __log2i(__last - __first)
    let depth = 2 * log2i(n) as isize;
    introsort(v, 0, n, &mut less, depth, true);
}

/// `std::partial_sort(first, first + middle, last, less)`.
pub fn partial_sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], middle: usize, mut less: F) {
    let n = v.len();
    if middle == 0 { return; }
    partial_sort_impl(v, 0, middle, n, &mut less);
}

/// Index of `std::max_element(first, last, less)`: the first maximum.
pub fn first_max_by<T, F: FnMut(&T, &T) -> bool>(v: &[T], mut less: F) -> Option<usize> {
    if v.is_empty() { return None; }
    let mut first = 0;
    for i in 1..v.len() {
        if less(&v[first], &v[i]) { first = i; }
    }
    Some(first)
}

/// Index of `std::min_element(first, last, less)` (also what `std::ranges::min` returns a copy of): the first minimum.
pub fn first_min_by<T, F: FnMut(&T, &T) -> bool>(v: &[T], mut less: F) -> Option<usize> {
    if v.is_empty() { return None; }
    let mut first = 0;
    for i in 1..v.len() {
        if less(&v[i], &v[first]) { first = i; }
    }
    Some(first)
}

/// `std::minmax(a, b)`: `(b, a)` when `b < a`, `(a, b)` otherwise.
#[inline]
pub fn minmax<T: PartialOrd>(a: T, b: T) -> (T, T) {
    if b < a { (b, a) } else { (a, b) }
}

/// `std::priority_queue<T, std::vector<T>, Compare>`: `top` is the maximum with respect to `less`
/// (`push` = push_back + push_heap, `pop` = pop_heap + pop_back).
pub struct PriorityQueue<T, F: FnMut(&T, &T) -> bool> {
    c: Vec<T>,
    comp: F,
}

impl<T, F: FnMut(&T, &T) -> bool> PriorityQueue<T, F> {
    pub fn new(less: F) -> Self { Self { c: Vec::new(), comp: less } }
    pub fn len(&self) -> usize { self.c.len() }
    pub fn is_empty(&self) -> bool { self.c.is_empty() }
    pub fn top(&self) -> Option<&T> { self.c.first() }

    pub fn push(&mut self, x: T) {
        self.c.push(x);
        let n = self.c.len();
        sift_up(&mut self.c, 0, n, &mut self.comp, n);
    }

    /// Removes and returns the top (C++ `top()` followed by `pop()`).
    pub fn pop(&mut self) -> Option<T> {
        let n = self.c.len();
        if n == 0 { return None; }
        pop_heap(&mut self.c, 0, n, &mut self.comp, n);
        self.c.pop()
    }
}

// ------------------------------------------------------------------------------------------------------------------

/// __log2i for a 64 bit difference_type: `sizeof(unsigned long) * CHAR_BIT - 1 - __libcpp_clz(n)`, 0 for n = 0.
fn log2i(n: usize) -> usize {
    if n == 0 { return 0; }
    63 - (n as u64).leading_zeros() as usize
}

// __sort3, non-branchless: stable, 2-3 compares, 0-2 swaps
fn sort3<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x: usize, y: usize, z: usize, c: &mut F) -> bool {
    if !c(&v[y], &v[x]) {          // if x <= y
        if !c(&v[z], &v[y]) {      // if y <= z
            return false;          // x <= y && y <= z
        }
        v.swap(y, z);              // x <= z && y < z
        if c(&v[y], &v[x]) {       // if x > y
            v.swap(x, y);          // x < y && y <= z
        }
        return true;
    }
    if c(&v[z], &v[y]) {           // x > y, if y > z
        v.swap(x, z);              // x < y && y < z
        return true;
    }
    v.swap(x, y);                  // x > y && y <= z
    if c(&v[z], &v[y]) {           // if y > z
        v.swap(y, z);              // x <= y && y < z
    }
    true
}

// __sort4, non-branchless
fn sort4<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x1: usize, x2: usize, x3: usize, x4: usize, c: &mut F) {
    sort3(v, x1, x2, x3, c);
    if c(&v[x4], &v[x3]) {
        v.swap(x3, x4);
        if c(&v[x3], &v[x2]) {
            v.swap(x2, x3);
            if c(&v[x2], &v[x1]) {
                v.swap(x1, x2);
            }
        }
    }
}

// __sort5, non-branchless
fn sort5<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], x1: usize, x2: usize, x3: usize, x4: usize, x5: usize, c: &mut F) {
    sort4(v, x1, x2, x3, x4, c);
    if c(&v[x5], &v[x4]) {
        v.swap(x4, x5);
        if c(&v[x4], &v[x3]) {
            v.swap(x3, x4);
            if c(&v[x3], &v[x2]) {
                v.swap(x2, x3);
                if c(&v[x2], &v[x1]) {
                    v.swap(x1, x2);
                }
            }
        }
    }
}

// __insertion_sort (guarded). The inserted element `t` stays at `i` while the slots it passes are scanned, then the
// run [h, i] is rotated: the same comparisons and the same result as the hole-based libc++ loop.
fn insertion_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    if first == last { return; }
    for i in first + 1..last {
        let j = i - 1;
        if c(&v[i], &v[j]) {
            let mut h = j;
            while h != first && c(&v[i], &v[h - 1]) { h -= 1; }
            v[h..=i].rotate_right(1);
        }
    }
}

// __insertion_sort_unguarded: assumes the element at first - 1 is not greater than any element of the range.
fn insertion_sort_unguarded<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    if first == last { return; }
    for i in first + 1..last {
        let j = i - 1;
        if c(&v[i], &v[j]) {
            let mut h = j;
            while c(&v[i], &v[h - 1]) { h -= 1; } // no bounds check, as in libc++
            v[h..=i].rotate_right(1);
        }
    }
}

// __insertion_sort_incomplete: gives up after 8 insertions; returns whether the range ended sorted.
fn insertion_sort_incomplete<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) -> bool {
    match last - first {
        0 | 1 => return true,
        2 => {
            if c(&v[last - 1], &v[first]) { v.swap(first, last - 1); }
            return true;
        }
        3 => { sort3(v, first, first + 1, last - 1, c); return true; }
        4 => { sort4(v, first, first + 1, first + 2, last - 1, c); return true; }
        5 => { sort5(v, first, first + 1, first + 2, first + 3, last - 1, c); return true; }
        _ => {}
    }
    let mut j = first + 2;
    sort3(v, first, first + 1, j, c);
    const LIMIT: u32 = 8;
    let mut count = 0;
    let mut i = j + 1;
    while i != last {
        if c(&v[i], &v[j]) {
            let mut h = j;
            while h != first && c(&v[i], &v[h - 1]) { h -= 1; }
            v[h..=i].rotate_right(1);
            count += 1;
            if count == LIMIT { return i + 1 == last; }
        }
        j = i;
        i += 1;
    }
    true
}

// __partition_with_equals_on_right: v[first] is the pivot; equivalent elements go right of it. Returns the pivot
// position and whether the range was already partitioned. The pivot stays at `begin` during the scan (libc++ moves
// it out; the slot is never read there).
fn partition_with_equals_on_right<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) -> (usize, bool) {
    let begin = first;
    let mut first = first;
    let mut last = last;
    // find the first element greater or equal to the pivot (guarded by the median of three)
    loop {
        first += 1;
        if !c(&v[first], &v[begin]) { break; }
    }
    // find the last element less than the pivot
    if begin == first - 1 {
        while first < last {
            last -= 1;
            if c(&v[last], &v[begin]) { break; }
        }
    } else {
        loop {
            last -= 1;
            if c(&v[last], &v[begin]) { break; }
        }
    }
    let already_partitioned = first >= last;
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if !c(&v[first], &v[begin]) { break; }
        }
        loop {
            last -= 1;
            if c(&v[last], &v[begin]) { break; }
        }
    }
    let pivot_pos = first - 1;
    if begin != pivot_pos { v.swap(begin, pivot_pos); }
    (pivot_pos, already_partitioned)
}

// __partition_with_equals_on_left: equivalent elements go left of the pivot. Returns the position after the pivot.
fn partition_with_equals_on_left<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) -> usize {
    let begin = first;
    let mut first = first;
    let mut last = last;
    if c(&v[begin], &v[last - 1]) {
        // guarded
        loop {
            first += 1;
            if c(&v[begin], &v[first]) { break; }
        }
    } else {
        loop {
            first += 1;
            if first >= last || c(&v[begin], &v[first]) { break; }
        }
    }
    if first < last {
        // guarded by the median of three
        loop {
            last -= 1;
            if !c(&v[begin], &v[last]) { break; }
        }
    }
    while first < last {
        v.swap(first, last);
        loop {
            first += 1;
            if c(&v[begin], &v[first]) { break; }
        }
        loop {
            last -= 1;
            if !c(&v[begin], &v[last]) { break; }
        }
    }
    let pivot_pos = first - 1;
    if begin != pivot_pos { v.swap(begin, pivot_pos); }
    first
}

// __introsort<_AlgPolicy, _Compare, _RandomAccessIterator, /*_UseBitSetPartition*/ false>
fn introsort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut first: usize, mut last: usize, c: &mut F, mut depth: isize, mut leftmost: bool) {
    // upper bound for using insertion sort
    const LIMIT: usize = 24;
    // lower bound for using Tukey's ninther for the pivot
    const NINTHER_THRESHOLD: usize = 128;
    loop {
        let len = last - first;
        match len {
            0 | 1 => return,
            2 => {
                if c(&v[last - 1], &v[first]) { v.swap(first, last - 1); }
                return;
            }
            3 => { sort3(v, first, first + 1, last - 1, c); return; }
            4 => { sort4(v, first, first + 1, first + 2, last - 1, c); return; }
            5 => { sort5(v, first, first + 1, first + 2, first + 3, last - 1, c); return; }
            _ => {}
        }
        if len < LIMIT {
            if leftmost { insertion_sort(v, first, last, c); } else { insertion_sort_unguarded(v, first, last, c); }
            return;
        }
        if depth == 0 {
            // fallback to heap sort, as introsort suggests
            partial_sort_impl(v, first, last, last, c);
            return;
        }
        depth -= 1;
        {
            let half = len / 2;
            if len > NINTHER_THRESHOLD {
                sort3(v, first, first + half, last - 1, c);
                sort3(v, first + 1, first + (half - 1), last - 2, c);
                sort3(v, first + 2, first + (half + 1), last - 3, c);
                sort3(v, first + (half - 1), first + half, first + (half + 1), c);
                v.swap(first, first + half);
            } else {
                sort3(v, first + half, first, last - 1, c);
            }
        }
        // the elements left of the range are sorted; if the pivot equals the greatest of them, every element of
        // [first, pivot] equals the pivot and the left side needs no sorting
        if !leftmost && !c(&v[first - 1], &v[first]) {
            first = partition_with_equals_on_left(v, first, last, c);
            continue;
        }
        let (i, already_partitioned) = partition_with_equals_on_right(v, first, last, c);
        // given a perfect partition, see whether insertion sort finishes quickly
        if already_partitioned {
            let fs = insertion_sort_incomplete(v, first, i, c);
            if insertion_sort_incomplete(v, i + 1, last, c) {
                if fs { return; }
                last = i;
                continue;
            } else if fs {
                first = i + 1;
                continue;
            }
        }
        // the left partition recursively, the right one by tail recursion elimination
        introsort(v, first, i, c, depth, leftmost);
        leftmost = false;
        first = i + 1;
    }
}

// __partial_sort_impl(first, middle, last)
fn partial_sort_impl<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, middle: usize, last: usize, c: &mut F) {
    if first == middle { return; }
    make_heap(v, first, middle, c);
    let len = middle - first;
    for i in middle..last {
        if c(&v[i], &v[first]) {
            v.swap(i, first);
            sift_down(v, first, len, first, c);
        }
    }
    sort_heap(v, first, middle, c);
}

// __make_heap
fn make_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    let n = last - first;
    if n > 1 {
        let mut start = (n - 2) / 2;
        loop {
            sift_down(v, first, n, first + start, c);
            if start == 0 { break; }
            start -= 1;
        }
    }
}

// __sort_heap
fn sort_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F) {
    let mut last = last;
    let mut n = last - first;
    while n > 1 {
        pop_heap(v, first, last, c, n);
        last -= 1;
        n -= 1;
    }
}

// __sift_down: the value at `start` (libc++'s `__top`) travels down by swaps.
fn sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, len: usize, start: usize, c: &mut F) {
    let mut child = start - first;
    if len < 2 || (len - 2) / 2 < child { return; }
    child = 2 * child + 1;
    let mut child_i = first + child;
    if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
        child_i += 1;
        child += 1;
    }
    if c(&v[child_i], &v[start]) { return; }
    let mut start = start; // `__top` sits here
    loop {
        v.swap(start, child_i);
        start = child_i;
        if (len - 2) / 2 < child { break; }
        child = 2 * child + 1;
        child_i = first + child;
        if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        if c(&v[child_i], &v[start]) { break; } // while (!__comp(*__child_i, __top))
    }
}

// __floyd_sift_down: the top travels to a leaf; returns its position (the hole).
fn floyd_sift_down<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, len: usize, c: &mut F) -> usize {
    debug_assert!(len >= 2);
    let mut hole = first;
    let mut child_i = first;
    let mut child = 0usize;
    loop {
        child_i += child + 1;
        child = 2 * child + 1;
        if child + 1 < len && c(&v[child_i], &v[child_i + 1]) {
            child_i += 1;
            child += 1;
        }
        v.swap(hole, child_i);
        hole = child_i;
        if child > (len - 2) / 2 { return hole; }
    }
}

// __pop_heap
fn pop_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F, len: usize) {
    if len > 1 {
        let hole = floyd_sift_down(v, first, len, c); // the top now sits at `hole`
        let last = last - 1;
        if hole != last {
            v.swap(hole, last); // *hole = move(*last); *last = move(top)
            let hole = hole + 1;
            sift_up(v, first, hole, c, hole - first);
        }
    }
}

// __sift_up: the value at last - 1 (libc++'s `__t`) travels up by swaps.
fn sift_up<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, c: &mut F, len: usize) {
    if len > 1 {
        let mut len = (len - 2) / 2;
        let mut ptr = first + len;
        let mut last = last - 1;
        if c(&v[ptr], &v[last]) {
            loop {
                v.swap(last, ptr);
                last = ptr;
                if len == 0 { break; }
                len = (len - 1) / 2;
                ptr = first + len;
                if !c(&v[ptr], &v[last]) { break; } // while (__comp(*__ptr, __t))
            }
        }
    }
}

// Expected outputs below were produced by the reference libc++ (Apple clang 17, -std=c++23 -O2) on the same inputs.
#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Copy)]
    struct S { key: i32, id: i32 }

    #[test]
    fn sort_lambda_with_ties() { // std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a.key < b.key; })
        let mut v: Vec<S> = (0..40).map(|i| S { key: (i * 7 + 3) % 4, id: i }).collect();
        sort_by(&mut v, |a, b| a.key < b.key);
        let ids: Vec<i32> = v.iter().map(|s| s.id).collect();
        assert_eq!(ids, [11, 31, 27, 3, 23, 39, 19, 7, 15, 35, 30, 10, 14, 34, 18, 6, 22, 38, 26, 2, 9, 13, 37, 17, 5, 21, 33, 25, 1, 29, 20, 28, 32, 24, 16, 12, 36, 8, 4, 0]);
    }

    #[test]
    fn ranges_sort_projection_ninther() { // rg::sort(v, {}, [&](int i) { return key[i]; }) on a vector<int>
        let key: Vec<i32> = (0..200).map(|i: i32| (i * i * 31 + 17 * i) % 5).collect();
        let mut v: Vec<i32> = (0..200).collect();
        sort_by(&mut v, |&a, &b| key[a as usize] < key[b as usize]);
        assert_eq!(v, [80, 198, 195, 3, 193, 5, 190, 188, 8, 185, 10, 183, 180, 13, 178, 15, 175, 173, 18, 170, 20, 168, 165, 23, 163, 25, 160, 158, 28, 155, 30, 153, 150, 33, 148, 35, 145, 143, 38, 140, 40, 138, 135, 43, 133, 45, 130, 128, 48, 125, 50, 123, 120, 53, 118, 55, 115, 113, 58, 110, 60, 108, 105, 63, 103, 65, 0, 100, 68, 98, 70, 95, 93, 73, 90, 75, 88, 85, 78, 83, 1, 161, 82, 197, 196, 77, 86, 87, 76, 2, 192, 91, 92, 72, 191, 71, 96, 97, 6, 67, 66, 101, 102, 7, 187, 62, 106, 107, 61, 186, 11, 111, 112, 57, 182, 56, 116, 117, 181, 12, 52, 121, 122, 51, 177, 176, 126, 127, 47, 16, 46, 131, 132, 17, 172, 42, 136, 137, 41, 171, 81, 141, 142, 37, 21, 36, 146, 147, 167, 166, 32, 151, 152, 31, 22, 162, 156, 157, 27, 26, 39, 199, 29, 24, 164, 154, 149, 34, 144, 169, 19, 139, 134, 44, 174, 129, 49, 124, 14, 179, 119, 54, 114, 59, 184, 9, 109, 104, 64, 189, 69, 94, 74, 4, 194, 89, 84, 79, 99, 159]);
    }

    #[test]
    fn sort_inconsistent_comparator() { // the rank-then-cross comparator of motorcycle.h (not a strict weak order)
        let rank = |i: i32| (i * 5) % 2;
        let dx = |i: i32| (i * 3) % 7 - 3;
        let dy = |i: i32| (i * 5) % 7 - 3;
        let mut v: Vec<[i32; 2]> = (0..12).map(|i| [i, 0]).collect();
        sort_by(&mut v, |a, b| {
            let (ra, rb) = (rank(a[0]), rank(b[0]));
            if ra != rb { return ra < rb; }
            dx(a[0]) * dy(b[0]) - dy(a[0]) * dx(b[0]) > 0
        });
        let ids: Vec<i32> = v.iter().map(|p| p[0]).collect();
        assert_eq!(ids, [0, 6, 2, 4, 8, 10, 1, 5, 7, 3, 9, 11]);
    }

    #[test]
    fn sort_adversary_hits_heap_fallback() { // McIlroy's killer adversary: same comparison count and output as libc++
        let n = 300usize;
        let gas = n as i64;
        let mut val = vec![gas; n];
        let (mut nsolid, mut candidate, mut ncmp) = (0i64, 0usize, 0u64);
        let mut v: Vec<i32> = (0..n as i32).collect();
        sort_by(&mut v, |&x, &y| {
            ncmp += 1;
            let (x, y) = (x as usize, y as usize);
            if val[x] == gas && val[y] == gas {
                if x == candidate { val[x] = nsolid; } else { val[y] = nsolid; }
                nsolid += 1;
            }
            if val[x] == gas { candidate = x; } else if val[y] == gas { candidate = y; }
            val[x] < val[y]
        });
        let h = v.iter().fold(1469598103934665603u64, |h, &x| (h ^ x as u64).wrapping_mul(1099511628211));
        assert_eq!((ncmp, h), (5983, 12636784802174111653));
    }

    #[test]
    fn priority_queue_with_ties() { // std::priority_queue<S, vector<S>, decltype([](a, b) { return a.key > b.key; })>
        let mut q = PriorityQueue::new(|a: &S, b: &S| a.key > b.key);
        let (mut out, mut next) = (vec![], 0);
        for i in 0..60 {
            if i % 4 == 3 { out.push(q.pop().unwrap().id); } else { q.push(S { key: (i * 13) % 5, id: next }); next += 1; }
        }
        while let Some(s) = q.pop() { out.push(s.id); }
        assert_eq!(out, [0, 4, 8, 2, 9, 15, 19, 23, 13, 24, 30, 34, 38, 32, 28, 17, 39, 43, 26, 11, 22, 33, 3, 18, 37, 41, 7, 1, 12, 35, 27, 5, 31, 20, 16, 42, 10, 29, 44, 21, 25, 14, 6, 36, 40]);
    }

    #[test]
    fn partial_sort_with_ties() { // std::partial_sort(v.begin(), v.begin() + 10, v.end(), less on key)
        let mut v: Vec<S> = (0..50).map(|i| S { key: (i * 11) % 7, id: i }).collect();
        partial_sort_by(&mut v, 10, |a, b| a.key < b.key);
        let ids: Vec<i32> = v.iter().map(|s| s.id).collect();
        assert_eq!(ids, [28, 21, 42, 0, 14, 49, 35, 7, 2, 9, 5, 3, 12, 10, 8, 15, 1, 17, 13, 19, 20, 6, 22, 11, 24, 25, 26, 27, 4, 29, 18, 31, 32, 33, 34, 16, 36, 37, 38, 39, 40, 41, 23, 43, 44, 45, 46, 47, 48, 30]);
    }

    #[test]
    fn first_extrema() { // std::max_element / std::min_element return the first of equal extrema
        let v = [1, 3, 3, 0, 0];
        assert_eq!(first_max_by(&v, |a, b| a < b), Some(1));
        assert_eq!(first_min_by(&v, |a, b| a < b), Some(3));
        assert_eq!(first_max_by::<i32, _>(&[], |a, b| a < b), None);
        assert_eq!(minmax(2, 1), (1, 2));
    }
}

/// `std::unordered_map<int, V>` of libc++ 20.1 (std::hash<int>, max_load_factor 1), reproducing its iteration order:
/// the C++ walks `per_face` below, and the order decides which contact of a pair of tedges is reported (its face and
/// the roles of the two tedges), hence which tedge `repair_crossing_tedges` reroutes first. Port of the insertion
/// (`__emplace_unique_key_args`) and rehash (`__rehash<true>`, `__do_rehash<true>`) of <__hash_table>: a singly
/// linked list of nodes, every bucket pointing at the node before its first one.
pub struct LibcxxIntMap<V> {
    keys: Vec<i32>,
    vals: Vec<V>,
    next: Vec<usize>,
    first: usize,        // __first_node_.__next_
    buckets: Vec<usize>, // the node before the bucket's first node (BEFORE_BEGIN: __first_node_), or NIL
    size: usize,
}

const NIL: usize = usize::MAX;
const BEFORE_BEGIN: usize = usize::MAX - 1;

fn libcxx_hash(k: i32) -> usize { k as i64 as usize } // std::hash<int>: static_cast<size_t>
fn is_hash_power2(bc: usize) -> bool { bc > 2 && bc & (bc - 1) == 0 }
fn constrain_hash(h: usize, bc: usize) -> usize {
    if bc & (bc - 1) == 0 { h & (bc - 1) } else if h < bc { h } else { h % bc }
}
fn next_hash_pow2(n: usize) -> usize { if n < 2 { n } else { 1usize << (usize::BITS - (n - 1).leading_zeros()) } }
/// `std::__next_prime(n)`: the smallest prime >= n
fn next_prime(n: usize) -> usize {
    let is_prime = |m: usize| {
        if m < 2 { return false; }
        let mut d = 2;
        while d * d <= m { if m % d == 0 { return false; } d += 1; }
        true
    };
    let mut m = n;
    while !is_prime(m) { m += 1; }
    m
}

impl<V: Default> Default for LibcxxIntMap<V> {
    fn default() -> Self { Self::new() }
}

impl<V: Default> LibcxxIntMap<V> {
    pub fn new() -> Self { LibcxxIntMap { keys: vec![], vals: vec![], next: vec![], first: NIL, buckets: vec![], size: 0 } }

    fn next_of(&self, p: usize) -> usize { if p == BEFORE_BEGIN { self.first } else { self.next[p] } }
    fn set_next(&mut self, p: usize, v: usize) { if p == BEFORE_BEGIN { self.first = v; } else { self.next[p] = v; } }

    /// `operator[]`
    pub fn entry(&mut self, k: i32) -> &mut V {
        let h = libcxx_hash(k);
        let mut bc = self.buckets.len();
        if bc != 0 {
            let chash = constrain_hash(h, bc);
            let nd = self.buckets[chash];
            if nd != NIL {
                let mut nd = self.next_of(nd);
                while nd != NIL && (libcxx_hash(self.keys[nd]) == h || constrain_hash(libcxx_hash(self.keys[nd]), bc) == chash) {
                    if libcxx_hash(self.keys[nd]) == h && self.keys[nd] == k { return &mut self.vals[nd]; }
                    nd = self.next[nd];
                }
            }
        }
        let n = self.keys.len();
        self.keys.push(k);
        self.vals.push(V::default());
        self.next.push(NIL);
        if (self.size + 1) as f32 > bc as f32 * 1.0f32 || bc == 0 {
            let a = 2 * bc + usize::from(!is_hash_power2(bc));
            let b = ((self.size + 1) as f32 / 1.0f32).ceil() as usize;
            self.rehash_unique(a.max(b));
            bc = self.buckets.len();
        }
        let chash = constrain_hash(h, bc);
        let pn = self.buckets[chash];
        if pn == NIL {
            let pn = BEFORE_BEGIN;
            self.next[n] = self.next_of(pn);
            self.set_next(pn, n);
            self.buckets[chash] = pn;
            if self.next[n] != NIL { let c = constrain_hash(libcxx_hash(self.keys[self.next[n]]), bc); self.buckets[c] = n; }
        } else {
            self.next[n] = self.next_of(pn);
            self.set_next(pn, n);
        }
        self.size += 1;
        &mut self.vals[n]
    }

    fn rehash_unique(&mut self, mut n: usize) {
        if n == 1 { n = 2; } else if n & (n - 1) != 0 { n = next_prime(n); }
        let bc = self.buckets.len();
        if n > bc {
            self.do_rehash(n);
        } else if n < bc {
            let m = ((self.size as f32) / 1.0f32).ceil() as usize;
            n = n.max(if is_hash_power2(bc) { next_hash_pow2(m) } else { next_prime(m) });
            if n < bc { self.do_rehash(n); }
        }
    }

    fn do_rehash(&mut self, nbc: usize) {
        self.buckets = vec![NIL; nbc];
        if nbc == 0 { return; }
        let mut pp = BEFORE_BEGIN;
        let mut cp = self.next_of(pp);
        if cp == NIL { return; }
        let mut chash = constrain_hash(libcxx_hash(self.keys[cp]), nbc);
        self.buckets[chash] = pp;
        let mut phash = chash;
        pp = cp;
        cp = self.next[cp];
        while cp != NIL {
            chash = constrain_hash(libcxx_hash(self.keys[cp]), nbc);
            if chash == phash {
                pp = cp;
            } else if self.buckets[chash] == NIL {
                self.buckets[chash] = pp;
                pp = cp;
                phash = chash;
            } else {
                let np = cp;
                let after = self.next[np];
                self.set_next(pp, after);
                let b = self.buckets[chash];
                self.next[np] = self.next_of(b);
                self.set_next(b, cp);
            }
            cp = self.next_of(pp);
        }
    }

    /// iteration order of the map: `(key, value)` from `begin()` to `end()`
    pub fn iter(&self) -> impl Iterator<Item = (i32, &V)> + '_ {
        let mut p = self.first;
        std::iter::from_fn(move || {
            if p == NIL { return None; }
            let r = (self.keys[p], &self.vals[p]);
            p = self.next[p];
            Some(r)
        })
    }
}


#[cfg(test)]
mod libcxx_int_map_tests {
    use super::LibcxxIntMap;

    // iteration orders produced by libc++ 20.1 (Apple clang 17, SDK) for `std::unordered_map<int, int> m; for (k :
    // keys) m[k]++;` then `for (auto& [k, v] : m)`
    #[test]
    fn libcxx_int_map_order() {
        let cases: [(&[i32], &[i32]); 3] = [
            (&[5, 3, 9, 1, 3, 7, 12, 0, 5, 44], &[44, 0, 7, 12, 1, 9, 3, 5]),
            (&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], &[11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]),
            (&[100, 23, 46, 69, 92, 1, 24, 47], &[92, 47, 69, 24, 46, 1, 23, 100]),
        ];
        for (keys, expected) in cases {
            let mut m: LibcxxIntMap<i32> = LibcxxIntMap::new();
            for &k in keys { *m.entry(k) += 1; }
            let got: Vec<i32> = m.iter().map(|(k, _)| k).collect();
            assert_eq!(got, expected, "keys {keys:?}");
        }
    }
}
