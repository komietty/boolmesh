use crate::common::Row2f;
use crate::triangulation::Pt;

pub fn compute_flat_tree(pts: &mut [Pt]) {
    if pts.len() <= 8 { return; }
    compute_flat_tree_impl(pts, true);
}

fn compute_flat_tree_impl(pts: &mut [Pt], sort_x: bool) {
    let eq = std::cmp::Ordering::Equal;
    if sort_x { pts.sort_by(|a, b| a.pos.x.partial_cmp(&b.pos.x).unwrap_or(eq)); }
    else      { pts.sort_by(|a, b| a.pos.y.partial_cmp(&b.pos.y).unwrap_or(eq)); }

    if pts.len() < 2 { return; }

    let (l, mr) = pts.split_at_mut(pts.len() / 2);
    if !mr.is_empty() {
        let (_, r)  = mr.split_first_mut().unwrap();
        compute_flat_tree_impl(l, !sort_x);
        compute_flat_tree_impl(r, !sort_x);
    }
}

pub fn compute_query_flat_tree<F>(
    pts: &[Pt],
    rect: &Rect,
    mut func: F,
) where F: FnMut(&Pt) {
    //for p in pts.iter() { if rect.contains(&p.pos) { func(p);} }
    //return;
    // Below is the more efficient implementation
    if pts.len() <= 8 {
        for p in pts.iter() { if rect.contains(&p.pos) { func(p);} }
    } else {
        query_two_d_tree(pts, rect.clone(), func);
    }
}

pub fn query_two_d_tree<F>(pts: &[Pt], r: Rect, mut f: F) where F: FnMut(&Pt) {
    let mut cur: Rect = Rect::default();
    let mut lev: i32 = 0;
    let mut bgn: usize = 0;
    let mut len: usize = pts.len();

    cur.min = Row2f::new(f64::MIN, f64::MIN);
    cur.max = Row2f::new(f64::MAX, f64::MAX);

    // Stack holds deferred right subtrees: (rect, start, len, level)
    let mut stack: Vec<(Rect, usize, usize, i32)> = Vec::with_capacity(64);

    loop {
        if len <= 8 {
            for i in 0..len {
                let p = &pts[bgn + i];
                if r.contains(&p.pos) { f(p); }
            }
            if let Some((r, b, ln, lv)) = stack.pop() {
                cur = r;
                bgn = b;
                len = ln;
                lev = lv;
                continue;
            } else { break; }
        }

        let mid_oft = len / 2;
        let mid_idx = bgn + mid_oft;
        let mid     = &pts[mid_idx];

        let mut rect_l = cur.clone();
        let mut rect_r = cur.clone();
        if lev % 2 == 0 {
            rect_l.max.x = mid.pos.x;
            rect_r.min.x = mid.pos.x;
        } else {
            rect_l.max.y = mid.pos.y;
            rect_r.min.y = mid.pos.y;
        }

        if r.contains(&mid.pos) { f(mid); }

        let overlaps_l = rect_l.overlap(&r);
        let overlaps_r = rect_r.overlap(&r);

        if overlaps_l {
            if overlaps_r {
                let r_bgn = mid_idx + 1;
                let r_len = len - (mid_oft + 1);
                stack.push((rect_r, r_bgn, r_len, lev + 1));
            }
            cur = rect_l;
            len = mid_oft;
            lev += 1;
        } else {
            cur = rect_r;
            bgn = mid_idx + 1;
            len = len - (mid_oft + 1);
            lev += 1;
        }
    }
}

#[derive(Clone)]
pub struct Rect {
    pub min: Row2f,
    pub max: Row2f,
}

impl Rect {
    pub fn default() -> Self {
        Self {
            min: Row2f::new(f64::MAX, f64::MAX),
            max: Row2f::new(f64::MIN, f64::MIN),
        }
    }
    pub fn new(a: &Row2f, b: &Row2f) -> Self {
        Self {
            min: Row2f::new(a.x.min(b.x), a.y.min(b.y)),
            max: Row2f::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    pub fn contains(&self, p: &Row2f) -> bool {
        p.x >= self.min.x &&
        p.x <= self.max.x &&
        p.y >= self.min.y &&
        p.y <= self.max.y
    }

    pub fn overlap(&self, r: &Rect) -> bool {
        self.max.x >= r.min.x &&
        self.max.y >= r.min.y &&
        self.min.x <= r.max.x &&
        self.min.y <= r.max.y
    }

    pub fn size(&self) -> Row2f { self.max - self.min }

    pub fn scale(&self) -> f64 {
        let a_min = self.min.x.abs().max(self.min.y.abs());
        let a_max = self.max.x.abs().max(self.max.y.abs());
        a_min.max(a_max)
    }

    pub fn union(&mut self, p: &Row2f) {
        self.min = Row2f::new(self.min.x.min(p.x), self.min.y.min(p.y));
        self.max = Row2f::new(self.max.x.max(p.x), self.max.y.max(p.y));
    }

}
