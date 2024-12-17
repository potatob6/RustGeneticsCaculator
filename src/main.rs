use std::{borrow::Borrow, cell::RefCell, collections::{hash_map::Keys, HashMap, HashSet}, default, fmt::Display, fs::exists, hash::{BuildHasherDefault, Hash}, hint::assert_unchecked, io::{self, stdin, stdout, Write}, process::{exit, id}, rc::Rc, sync::LazyLock, time::SystemTime, u8, usize, vec};

use bigdecimal::{BigDecimal, FromPrimitive};
use fraction::Fraction;
use nohash::{BuildNoHashHasher, IsEnabled, NoHashHasher};

mod fraction;

type NoHashSet<T> = HashMap<T, (), nohash::BuildNoHashHasher<T>>;

macro_rules! nohashset {
    () => {
        HashMap::<GeneGroup, (), nohash::BuildNoHashHasher<GeneGroup>>::with_hasher(BuildHasherDefault::default())
    };
    ($n: expr) => {
        HashMap::<GeneGroup, (), nohash::BuildNoHashHasher<GeneGroup>>::with_capacity_and_hasher($n, BuildHasherDefault::default())
    }
}

#[repr(u8)]
#[derive(PartialEq)]
#[derive(Clone, Copy)]
#[derive(Debug)]
enum Gene {
    G, Y, X, W, H, PLACEHOLDER
}

const GENE_SCORE: [(Gene, u8); 6] = [
    (Gene::G, 5),
    (Gene::H, 5),
    (Gene::Y, 5),
    (Gene::W, 6),
    (Gene::X, 6),
    (Gene::PLACEHOLDER, 0)];

#[inline(always)]
fn get_score(pattern: &Gene) -> (usize, u8) {
    if GENE_SCORE[0].0 == *pattern {
        return (0, GENE_SCORE[0].1);
    }
    if GENE_SCORE[1].0 == *pattern {
        return (1, GENE_SCORE[1].1);
    }
    if GENE_SCORE[2].0 == *pattern {
        return (2, GENE_SCORE[2].1);
    }
    if GENE_SCORE[3].0 == *pattern {
        return (3, GENE_SCORE[3].1);
    }
    if GENE_SCORE[4].0 == *pattern {
        return (4, GENE_SCORE[4].1);
    }
    panic!()
}

#[derive(Debug)]
struct GeneGroup([Gene; 6], usize);

impl Display for GeneGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for i in 0..self.0.len() {
            if self.0[i] == Gene::G {
                s += &"G";
                continue;
            }
            if self.0[i] == Gene::Y {
                s += &"Y";
                continue;
            }
            if self.0[i] == Gene::X {
                s += &"X";
                continue;
            }
            if self.0[i] == Gene::H {
                s += &"H";
                continue;
            }
            if self.0[i] == Gene::W {
                s += &"W";
                continue;
            }
        }
        f.write_str(&s).expect("Error");
        Ok(())
    }
}

#[derive(Debug)]
#[derive(Clone)]
struct ComposeResult {
    sign: usize,
    loss: RefCell<u8>,
    gene_group: GeneGroup,
    probability: Fraction,
    prev_gene_group: Vec<GeneGroup>,
}

impl Hash for GeneGroup {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Only use lower 48 bits
        state.write_u64(self.hash_expr());
    }
}

impl Hash for ComposeResult {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.gene_group.hash(state);
    }
}

impl Eq for ComposeResult {

}

impl Eq for GeneGroup {
    
}

impl IsEnabled for ComposeResult {

}

impl IsEnabled for GeneGroup {

}

impl PartialEq for ComposeResult {
    fn eq(&self, other: &Self) -> bool {
        let g1 = self.gene_group == other.gene_group && self.probability.1 == other.probability.1;
        if !g1 {
            return false;
        }
        
        if self.prev_gene_group.len() != other.prev_gene_group.len() {
            return false;
        }

        let result = true;
        for i in 0..self.prev_gene_group.len() {
            if self.prev_gene_group[i].1 != other.prev_gene_group[i].1 {
                return false;
            }

            if self.prev_gene_group[i].1 == 0 && other.prev_gene_group[i].1 == 0 {
                if self.prev_gene_group[i].0 != other.prev_gene_group[i].0 {
                    return false;
                }
            }
        }
        result
    }
}

impl PartialEq for GeneGroup {
    fn eq(&self, other: &Self) -> bool {
        if self.1 != 0 && other.1 != 0 && self.1 == other.1 {
            return true;
        }
        let mut result = true;
        for i in 0..6 {
            if self.0[i] != other.0[i] {
                result = false;
                break;
            }
        }
        result
    }
}

impl GeneGroup {

    #[inline(always)]
    fn hash_expr(&self) -> u64 {
        let result: u64 
        = ((self.0[0] as u64) << 40 | 
          (self.0[1] as u64) << 32 | 
          (self.0[2] as u64) << 24 | 
          (self.0[3] as u64) << 16 |
          (self.0[4] as u64) << 8 |
          (self.0[5] as u64));
        result
    }

    fn from_str(s: &str) -> Option<Self> {
        if s.len() < 6 {
            return None;
        }
        let mut s_result = Self([Gene::PLACEHOLDER; 6], 0);
        for (i, c) in s.chars().enumerate() {
            if c == 'X' || c == 'x' {
                s_result.0[i] = Gene::X;
            } else if c == 'Y' || c == 'y' {
                s_result.0[i] = Gene::Y;
            } else if c == 'H' || c == 'h' {
                s_result.0[i] = Gene::H;
            } else if c == 'W' || c == 'w' {
                s_result.0[i] = Gene::W;
            } else if c == 'G' || c == 'g' {
                s_result.0[i] = Gene::G;
            } else {
                return None;
            }
        }

        Some(s_result)
    }

    #[inline]
    fn display(&self) -> String {
        let mut s = String::new();
        for i in 0..6 {
            if self.0[i] == Gene::X {
                s.push('X');
                continue;
            }
            if self.0[i] == Gene::Y {
                s.push('Y');
                continue;
            }
            if self.0[i] == Gene::G {
                s.push('G');
                continue;
            }
            if self.0[i] == Gene::H {
                s.push('H');
                continue;
            }
            if self.0[i] == Gene::W {
                s.push('W');
                continue;
            }
        }
        s
    }
}

impl Clone for GeneGroup {
    fn clone(&self) -> Self {
        Self(self.0, self.1)
    }
}

static mut GLOBAL_SIGN: usize = 1;

fn compose(plants: Vec<&GeneGroup>, exists_gene: &NoHashSet<GeneGroup>, 
    already: &NoHashSet<ComposeResult>, selected: &NoHashSet<ComposeResult>, 
    probability: &BigDecimal) -> Vec<ComposeResult> {

    unsafe { assert_unchecked(plants.len() >= 2 && plants.len() <= 4); }
    let mut all_situation: [Vec<Gene>; 6] = [const { Vec::new() }; 6];
    for i in 0..6 {
        let mut scores = [0u8; 6];
        for plant in &plants {
            let (id, s) = get_score(&plant.0[i]);
            scores[id] += s;
        }

        // Get max gene score
        let mut max_score = 0u8;
        for o in 0..6 {
            if scores[o] > max_score {
                max_score = scores[o];
            }
        }

        // Get all max score gene
        let mut max_gene_score = Vec::<Gene>::new();
        for o in 0..6 {
            if scores[o] == max_score {
                max_gene_score.push(GENE_SCORE[o].0);
            }
        }

        all_situation[i] = max_gene_score;
    }

    // Expand all situation
    let expanded = expand(&all_situation);

    expanded
        .into_iter()
        .filter_map(
            move |a,| { 
                compose_filter(a.0, a.1, probability, exists_gene, already, selected, &plants) 
            }
        )
        .collect::<Vec<ComposeResult>>()
}

#[inline(always)]
fn compose_filter(gene_group: GeneGroup, probability: Fraction, global_probability: &BigDecimal, 
                exists: &NoHashSet<GeneGroup>, already: &NoHashSet<ComposeResult>, selected: &NoHashSet<ComposeResult>, prev_gene_group: &Vec<&GeneGroup>) -> Option<ComposeResult> {
    if &probability.1 > global_probability {
        return None;
    }

    let g = ComposeResult {
        sign: 0,
        loss: RefCell::new(u8::MAX),
        gene_group: gene_group,
        probability: probability,
        prev_gene_group: prev_gene_group.iter().cloned().map(|x| x.clone() ).collect::<Vec<GeneGroup>>(),
    };

    if exists.contains_key(&g.gene_group) || already.contains_key(&g) || selected.contains_key(&g) {
        return None;
    }

    Some(g)
}

fn expand(sit: &[Vec<Gene>; 6]) -> Vec<(GeneGroup, Fraction)> {
    let mut expanded = Vec::<(GeneGroup, Fraction)>::new();
    #[derive(Debug)]
    struct StackNode {
        idx: usize,
        now_node: Gene,
        prev_gene_seq: Option<Vec<Gene>>,
        probability: Fraction,
    }

    let mut stack = Vec::<StackNode>::new();
    let mut leaf_path = Vec::<StackNode>::new();

    loop {
        if stack.len() == 0 {
            unsafe { assert_unchecked( stack.len() == 0 ); }
            for j in &sit[0] {
                let g = StackNode {
                    idx: 0,
                    now_node: j.clone(),
                    prev_gene_seq: None,
                    probability: Fraction::from_u32(1, (&sit[0]).len() as u32),
                };
                stack.push(g);
            }

            if stack.len() == 0 {
                break;
            }

            continue;
        }

        
        unsafe { assert_unchecked( stack.len() >= 1 ); }
        let stack_node = stack.pop().unwrap();
        let now_deal_idx = stack_node.idx + 1;
        let mut prev_seq = Vec::new();

        match &stack_node.prev_gene_seq {
            Some(_) => {
                unsafe { assert_unchecked(stack_node.prev_gene_seq.is_some()); }
                let mut prev_seq_1 = stack_node.prev_gene_seq.unwrap();
                prev_seq_1.push(stack_node.now_node);
                prev_seq = prev_seq_1;
            },
            None => {
                prev_seq = vec![stack_node.now_node];
            }
        }


        for j in &sit[now_deal_idx] {
            let g = StackNode {
                idx: stack_node.idx + 1,
                now_node: j.clone(),
                prev_gene_seq: Some(prev_seq.clone()),
                probability: &stack_node.probability * &Fraction::from_u32(1, (&sit[now_deal_idx]).len() as u32),
            };

            if now_deal_idx == 5 {
                leaf_path.push(g);
            } else {
                stack.push(g);
            }
        }

        if stack.len() == 0 {
            break;
        }
    }

    for leaf in leaf_path {
        unsafe { assert_unchecked(leaf.prev_gene_seq.is_some()); }
        let gene_seq = leaf.prev_gene_seq.unwrap();

        unsafe { assert_unchecked(gene_seq.len() == 5); }

        let g = GeneGroup([
            gene_seq[0], 
            gene_seq[1],
            gene_seq[2], 
            gene_seq[3],
            gene_seq[4],
            leaf.now_node
            ], 0);
            expanded.push((g, leaf.probability));
    }

    expanded
}

fn select2compose(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>) -> Vec<ComposeResult> {

    // already_collection only for check repeat item
    let mut output = Vec::with_capacity( exist_gene_vec.len().pow(2) );
    let l = exist_gene_vec.len();
    for i in 0..l {
        for j in i..l {
            if i != j {
                let v = vec![exist_gene_vec[i], exist_gene_vec[j]];
                let mut results = compose(v, exist_hash, already_collection, selected, probability);
                // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                output.append(&mut results);
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select3compose(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>) -> Vec<ComposeResult> {

    // already_collection only for check repeat item
    let mut output = Vec::with_capacity( exist_gene_vec.len().pow(3) );
    let l = exist_gene_vec.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                if !(i == j && i == k) {
                    let v = vec![exist_gene_vec[i], exist_gene_vec[j], exist_gene_vec[k]];
                    let mut results = compose(v, exist_hash, already_collection, selected, probability);
                    // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                    output.append(&mut results);
                }
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select4compose(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>) -> Vec<ComposeResult> {

    // already_collection only for check repeat item
    let mut output = Vec::with_capacity( exist_gene_vec.len().pow(4) );
    let l = exist_gene_vec.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                for z in k..l {
                    if !(i == j && i == k && i == z) {
                        let v = vec![exist_gene_vec[i], exist_gene_vec[j], exist_gene_vec[k], exist_gene_vec[z]];
                        let mut results = compose(v, exist_hash, already_collection, selected, probability);
                        // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                        output.append(&mut results);
                    }
                }
            }
        }
    }
    output
}

fn select2compose_delta(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>, fixed_node: &GeneGroup, 
    delta_compose: &mut NoHashSet<ComposeResult>) {

    // already_collection only for check repeat item
    // let mut output = Vec::with_capacity( exist_gene_vec.len().pow(2) );
    let l = exist_gene_vec.len();
    for i in 0..l {
        if exist_gene_vec[i] != fixed_node {
            let v = vec![exist_gene_vec[i], fixed_node];
            let results = compose(v, exist_hash, already_collection, selected, probability);
            // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
            for mut item in results {
                item.sign = unsafe { GLOBAL_SIGN };
                item.gene_group.1 = unsafe { GLOBAL_SIGN };
                unsafe { GLOBAL_SIGN += 1 };
                delta_compose.insert(item, ());
            }
        }
    }
}

/// Output is the added items for already_collection
fn select3compose_delta(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>, fixed_node: &GeneGroup, 
    delta_compose: &mut NoHashSet<ComposeResult>) {

    // already_collection only for check repeat item
    let l = exist_gene_vec.len();
    for i in 0..l {
        for j in i..l {
            if !(i == j && exist_gene_vec[i] == fixed_node) {
                let v = vec![exist_gene_vec[i], exist_gene_vec[j], fixed_node];
                let results = compose(v, exist_hash, already_collection, selected, probability);
                // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                for mut item in results {
                    item.sign = unsafe { GLOBAL_SIGN };
                    item.gene_group.1 = unsafe { GLOBAL_SIGN };
                    unsafe { GLOBAL_SIGN += 1 };
                    delta_compose.insert(item, ());
                }
            }
        }
    }
}

/// Output is the added items for already_collection
fn select4compose_delta(exist_gene_vec: &Vec<&GeneGroup>, exist_hash: &NoHashSet<GeneGroup>, 
    already_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    selected: &NoHashSet<ComposeResult>, fixed_node: &GeneGroup, 
    delta_compose: &mut NoHashSet<ComposeResult>) {

    // already_collection only for check repeat item
    let l = exist_gene_vec.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                if !(i == j && i == k && exist_gene_vec[i] == fixed_node) {
                    let v = vec![exist_gene_vec[i], exist_gene_vec[j], exist_gene_vec[k], fixed_node];
                    let results = compose(v, exist_hash, already_collection, selected, probability);
                    // let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                    for mut item in results {
                        item.sign = unsafe { GLOBAL_SIGN };
                        item.gene_group.1 = unsafe { GLOBAL_SIGN };
                        unsafe { GLOBAL_SIGN += 1 };
                        delta_compose.insert(item, ());
                    }
                }
            }
        }
    }
}

fn one_step_compose_predict(exists_gene: &NoHashSet<GeneGroup>, already_compose_collection: &NoHashSet<ComposeResult>, selected_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal) 
    -> (Vec<ComposeResult>, Vec<ComposeResult>, Vec<ComposeResult>) {

    let exists_vec = exists_gene.keys().collect::<Vec<&GeneGroup>>();

    let mut a1 = select2compose(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection);
    let mut a2 = select3compose(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection);
    let mut a3 = select4compose(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection);

    // let mut a2 = vec_result_non_exists(&a2, &a1, &already_compose_collection, probability);
    // a1.append(&mut a2);
    // let mut a3 = vec_result_non_exists(&a3, &a1, &already_compose_collection, probability);
    // a1.append(&mut a3);
    // a1 = vec_result_non_selected(a1, selected_collection, probability);

    for j in &mut a1 {
        // SAFTY: Use static data with single thread is safe. No any data conflition.
        j.sign = unsafe { GLOBAL_SIGN };
        j.gene_group.1 = unsafe { GLOBAL_SIGN };
        unsafe { GLOBAL_SIGN += 1 };
    }

    evaluate_loss(&mut a1);
    (a1, a2, a3)
}

fn one_step_compose_predict_delta(exists_gene: &NoHashSet<GeneGroup>, already_compose_collection: &NoHashSet<ComposeResult>, 
    selected_collection: &NoHashSet<ComposeResult>, probability: &BigDecimal, 
    fixed_node: &GeneGroup, delta_compose: &mut NoHashSet<ComposeResult>) {

    let exists_vec = exists_gene.keys().collect::<Vec<&GeneGroup>>();

    select2compose_delta(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection, fixed_node, delta_compose);
    select3compose_delta(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection, fixed_node, delta_compose);
    select4compose_delta(&exists_vec, exists_gene, already_compose_collection, probability, selected_collection, fixed_node, delta_compose);

    evaluate_loss_delta(delta_compose);
}

fn evaluate_loss<'a>(compose_result: &'a mut Vec<ComposeResult>) {
    let mut lowest_idx = 0;
    let mut lowest_loss_score: Option<&'a ComposeResult> = None;
    // SAFTY: Use static data with single thread is safe. No any data conflition.
    if unsafe { &LOSS_MODE } == &LossMode::ACCURACY {
        for idx in 0..compose_result.len() {
            let mut loss_score = 0u8;
            for j in 0..6 {
                if compose_result[idx].gene_group.0[j] != unsafe { ACCURACY_TARGET }[j] {
                    loss_score += 1;
                }
            }

            {
                *compose_result[idx].loss.borrow_mut() = loss_score;
            }

            match lowest_loss_score {
                Some(prev) => {
                    if loss_score < *prev.loss.borrow() {
                        lowest_loss_score = Some(&compose_result[idx]);
                        lowest_idx = idx;
                    } else {
                        lowest_loss_score = Some(prev);
                    }
                }
                None => {
                    lowest_loss_score = Some(&compose_result[idx]);
                    lowest_idx = idx;
                }
            }
        }

        return;
    } else if unsafe { &LOSS_MODE } == &LossMode::AMOUNT {
        let mut lowest_loss_score: Option<&'a ComposeResult> = None;
        for idx in 0..compose_result.len() {
            let mut seperate_gene_amount = [0u8; 5];
            let mut loss_score = 0u8;

            for k in 0..6 {

                let g = compose_result[idx].gene_group.0[k];

                // Count Gene::G
                if g == Gene::G {
                    seperate_gene_amount[0] += 1;
                    continue;
                }

                // Count Gene::H
                if g == Gene::H {
                    seperate_gene_amount[1] += 1;
                    continue;
                }

                // Count Gene::Y
                if g == Gene::Y {
                    seperate_gene_amount[2] += 1;
                    continue;
                }

                // Count Gene::W
                if g == Gene::W {
                    seperate_gene_amount[3] += 1;
                    continue;
                }

                // Count Gene::X
                if g == Gene::X {
                    seperate_gene_amount[4] += 1;
                    continue;
                }
            }

            for k in 0..5 {
                if seperate_gene_amount[k] >= unsafe { AMOUNT_TARGET }[k] {
                    loss_score += seperate_gene_amount[k] - unsafe { AMOUNT_TARGET }[k];
                } else {
                    loss_score += unsafe { AMOUNT_TARGET }[k] - seperate_gene_amount[k];
                }
            }

            { *compose_result[idx].loss.borrow_mut() = loss_score; }

            match lowest_loss_score {
                Some(prev) => {
                    if loss_score < *prev.loss.borrow() {
                        lowest_loss_score = Some(&compose_result[idx]);
                        lowest_idx = idx;
                    } else {
                        lowest_loss_score = Some(prev);
                    }
                }
                None => {
                    lowest_loss_score = Some(&compose_result[idx]);
                    lowest_idx = idx;
                }
            }
        }
        
        
        return;
    }
    panic!("Unknown Loss mode.");
}

fn evaluate_loss_delta<'a>(compose_result: &'a mut NoHashSet<ComposeResult>) -> Option<&'a ComposeResult> {
    let mut lowest_loss_score: Option<&'a ComposeResult> = None;
    if compose_result.len() == 0 {
        return None;
    }
    // SAFTY: Use static data with single thread is safe. No any data conflition.
    if unsafe { &LOSS_MODE } == &LossMode::ACCURACY {
        for (item, _) in compose_result {
            let mut loss_score = 0u8;
            for j in 0..6 {
                if item.gene_group.0[j] != unsafe { ACCURACY_TARGET }[j] {
                    loss_score += 1;
                }
            }

            {
                *item.loss.borrow_mut() = loss_score;
            }

            match lowest_loss_score {
                Some(prev) => {
                    if loss_score < *prev.loss.borrow() {
                        lowest_loss_score = Some(item);
                    } else {
                        lowest_loss_score = Some(prev);
                    }
                }
                None => {
                    lowest_loss_score = Some(item);
                }
            }
        }

        return lowest_loss_score;
    } else if unsafe { &LOSS_MODE } == &LossMode::AMOUNT {
        let mut lowest_loss_score: Option<&'a ComposeResult> = None;
        for (item, _) in compose_result {
            let mut seperate_gene_amount = [0u8; 5];
            let mut loss_score = 0u8;

            for k in 0..6 {

                let g = item.gene_group.0[k];

                // Count Gene::G
                if g == Gene::G {
                    seperate_gene_amount[0] += 1;
                    continue;
                }

                // Count Gene::H
                if g == Gene::H {
                    seperate_gene_amount[1] += 1;
                    continue;
                }

                // Count Gene::Y
                if g == Gene::Y {
                    seperate_gene_amount[2] += 1;
                    continue;
                }

                // Count Gene::W
                if g == Gene::W {
                    seperate_gene_amount[3] += 1;
                    continue;
                }

                // Count Gene::X
                if g == Gene::X {
                    seperate_gene_amount[4] += 1;
                    continue;
                }
            }

            for k in 0..5 {
                if seperate_gene_amount[k] >= unsafe { AMOUNT_TARGET }[k] {
                    loss_score += seperate_gene_amount[k] - unsafe { AMOUNT_TARGET }[k];
                } else {
                    loss_score += unsafe { AMOUNT_TARGET }[k] - seperate_gene_amount[k];
                }
            }

            { *item.loss.borrow_mut() = loss_score; }

            match lowest_loss_score {
                Some(prev) => {
                    if loss_score < *prev.loss.borrow() {
                        lowest_loss_score = Some(item);
                    } else {
                        lowest_loss_score = Some(prev);
                    }
                }
                None => {
                    lowest_loss_score = Some(item);
                }
            }
        }
        
        
        return lowest_loss_score;
    }
    panic!("Unknown Loss mode.");
}

fn first_time_step(exist_gene: &mut NoHashSet<GeneGroup>, composed_collection: &mut NoHashSet<ComposeResult>, selected_collection: &mut NoHashSet<ComposeResult>, probability: &BigDecimal)
    -> (Option<ComposeResult>, Option<ComposeResult>) {
    let (a1, a2, a3) = one_step_compose_predict(&exist_gene, &composed_collection, &selected_collection, probability);
    for i in a3 {
        composed_collection.insert(i, ());
    }
    for i in a2 {
        composed_collection.insert(i, ());
    }
    for i in a1 {
        composed_collection.insert(i, ());
    }

    if composed_collection.len() == 0 {
        return (None, None);
    }

    unsafe { assert_unchecked((&composed_collection).len() > 0); }
    // Find lowest loss
    let (mut lowest_loss, _) = composed_collection.iter().next().unwrap();
    for (k, _) in composed_collection.iter() {
        if *k.loss.borrow() < *lowest_loss.loss.borrow() {
            lowest_loss = k;
        }
    }

    let lowest_loss = lowest_loss.clone();

    if *lowest_loss.loss.borrow() != 0 {
        exist_gene.insert(lowest_loss.gene_group.clone(), ());
        // println!("{} lowest {}", exist_gene.len(), lowest_loss.clone().gene_group.display());
        selected_collection.insert(lowest_loss.clone(), ());
        composed_collection.remove(&lowest_loss);

        return (None, Some(lowest_loss.clone()));
    } else {
        return (Some(lowest_loss.clone()), Some(lowest_loss.clone()));
    }
}

fn delta_step(exist_gene: &mut NoHashSet<GeneGroup>, composed_collection: &mut NoHashSet<ComposeResult>, 
    selected_collection: &mut NoHashSet<ComposeResult>, probability: &BigDecimal, 
    last_selected_genegroup: &GeneGroup)
-> (Option<ComposeResult>, Option<ComposeResult>, NoHashSet<ComposeResult>) {

    let mut delta_compose = NoHashSet::<ComposeResult>::with_capacity_and_hasher(exist_gene.len().pow(4), BuildNoHashHasher::default());
    one_step_compose_predict_delta(&exist_gene, &composed_collection, &selected_collection, probability, last_selected_genegroup, &mut delta_compose);
    if delta_compose.len() == 0 {
        return (None, None, delta_compose);
    }

    unsafe { assert_unchecked((&composed_collection).len() > 0); }

    let (min, _) = delta_compose.iter()
        .min_by_key(|x| { *x.0.loss.borrow() }).unwrap();
    let min = min.clone();

    delta_compose.remove(&min.clone());

    if *min.loss.borrow() == 0 {
        return (Some(min), None, delta_compose);
    }

    return (None, Some(min), delta_compose);
}

fn predict(exist_gene: &mut NoHashSet<GeneGroup>, composed_collection: &mut NoHashSet<ComposeResult>, selected_collection: &mut NoHashSet<ComposeResult>, probability: &BigDecimal) -> (Option<ComposeResult>, Option<ComposeResult>) {
    let mut generation = 1usize;
    let mut partial_lowest: Option<ComposeResult> = None;
    
    let (global, mut partial) = first_time_step(exist_gene, composed_collection, selected_collection, probability);
    match global {
        Some(g) => {
            unsafe { assert_unchecked(partial.is_some()); }
            return (Some(g), partial);
        }
        None => { }
    }

    match partial {
        None => {
            return (None, None);
        }
        Some(mut prev_selected_result) => {
            // Delta update
            let mut prev_lowest_loss = composed_collection.iter().min_by_key(|x| *x.0.loss.borrow() );
            let mut prev_selected_genegroup = prev_selected_result.gene_group.clone();

            let mut prev_lowest_loss = match prev_lowest_loss {
                Some(p) => Some((p.0.clone(), ())),
                None => None,
            };
            loop {
                let (global, partial, mut delta_compose) = delta_step(exist_gene, composed_collection, selected_collection, probability, &prev_selected_genegroup);
                let delta_compose_lowest = delta_compose.iter().min_by_key(|x| *x.0.loss.borrow() );
                match global {
                    Some(_) => {
                        return (global, None);
                    }
                    None => { },
                }

                match partial {
                    None => { return (None, None) }
                    Some(now) => {
                        match prev_lowest_loss {
                            Some(ref prev_lowest_loss1) => {
                                if *now.loss.borrow() < *prev_lowest_loss1.0.loss.borrow() {
                                    prev_selected_genegroup = now.gene_group.clone();
                                    exist_gene.insert(now.gene_group.clone(), ());
                                    selected_collection.insert(now, ());

                                    match delta_compose_lowest {
                                        Some(delta_compose_lowest1) => {
                                            if *delta_compose_lowest1.0.loss.borrow() < *prev_lowest_loss1.0.loss.borrow() {
                                                prev_lowest_loss = Some((delta_compose_lowest1.0.clone(), ()));
                                            } else { }
                                        }
                                        None => { }
                                    }

                                    let _ = delta_compose.into_iter().map(|x| composed_collection.insert(x.0, ()) );
                                } else {
                                    composed_collection.remove(&prev_lowest_loss1.0);
                                    exist_gene.insert((&prev_lowest_loss1.0.gene_group).clone(), ());
                                    prev_selected_genegroup = (&prev_lowest_loss1.0.gene_group).clone();
                                    selected_collection.insert(prev_lowest_loss1.0.clone(), ());

                                    let prev_min = composed_collection.iter().min_by_key(|x| {
                                        *x.0.loss.borrow()
                                    });

                                    match prev_min {
                                        Some(min) => {
                                            if *now.loss.borrow() < *min.0.loss.borrow() {
                                                prev_lowest_loss = Some((now.clone(), ()));
                                            } else {
                                                prev_lowest_loss = None;
                                            }
                                        }
                                        None => {
                                            prev_lowest_loss = Some((now.clone(), ()));
                                        }
                                    }
                                    let _ = delta_compose.into_iter().map(|x| composed_collection.insert(x.0, ()) );
                                }
                            },
                            None => {
                                exist_gene.insert(now.gene_group.clone(), ());
                                prev_selected_genegroup = now.gene_group.clone();
                                selected_collection.insert(now, ());
                                match delta_compose_lowest {
                                    Some(delta_compose_lowest) => {
                                        prev_lowest_loss = Some((delta_compose_lowest.0.clone(), ()));
                                    }
                                    None => { }
                                }

                                let _ = delta_compose.into_iter().map(|x| composed_collection.insert(x.0.clone(), ()));
                            },
                        }
                    }
                }                

                if exist_gene.len() > unsafe { SPEARD_LIMIT } {
                    return (None, Some(prev_selected_result));
                }
            }
        }
    }
    generation += 1;
}

#[repr(u8)]
#[derive(PartialEq)]
enum LossMode {
    ACCURACY,
    AMOUNT,
}

static mut LOSS_MODE: LossMode = LossMode::AMOUNT;
static mut ACCURACY_TARGET: [Gene; 6] = [Gene::G; 6];

///
///     (Gene::G, 5),
///     (Gene::H, 5),
///     (Gene::Y, 5),
///     (Gene::W, 6),
///     (Gene::X, 6),
static mut AMOUNT_TARGET: [u8; 5] = [3, 0, 3, 0, 0];
static mut SPEARD_LIMIT: usize = 30;

fn display_target() -> String {
    let mut result = "".to_string();
    if unsafe { AMOUNT_TARGET }[0] != 0 {
        result += &format!("{}G", unsafe { AMOUNT_TARGET[0] });
    }
    if unsafe { AMOUNT_TARGET }[1] != 0 {
        result += &format!("{}H", unsafe { AMOUNT_TARGET[1] });
    }
    if unsafe { AMOUNT_TARGET }[2] != 0 {
        result += &format!("{}Y", unsafe { AMOUNT_TARGET[2] });
    }
    if unsafe { AMOUNT_TARGET }[3] != 0 {
        result += &format!("{}W", unsafe { AMOUNT_TARGET[3] });
    }
    if unsafe { AMOUNT_TARGET }[4] != 0 {
        result += &format!("{}X", unsafe { AMOUNT_TARGET[4] });
    }
    result
}

fn find_prev_compose_node(node: ComposeResult, selected_collection: &NoHashSet<ComposeResult>) -> HashMap<usize, (ComposeResult, usize), BuildNoHashHasher<usize>> {
    let mut v = HashMap::<usize, (ComposeResult, usize), BuildNoHashHasher<usize>>::with_hasher(BuildNoHashHasher::default());
    for k in node.prev_gene_group {
        if k.1 != 0 {
            if !v.contains_key(&k.1) {
                // Find compose result
                for c in selected_collection {
                    if c.0.sign == k.1 {
                        v.insert(k.1, (c.0.clone(), 1));
                    }
                }
            } else {
                unsafe { assert_unchecked(v.get(&k.1).is_some()); }
                let g = v.get_mut(&k.1).unwrap();
                g.1 += 1;
            }
        }
    }
    v
}

fn display_backtrace_path(result: ComposeResult, selected_compose: &NoHashSet<ComposeResult>) -> String {

    #[derive(Debug)]
    struct BackTraceNode {
        times: usize, 
        node: ComposeResult,
    }


    let mut leaf = Vec::<(GeneGroup, usize)>::new();
    let mut output = String::new();
    let mut tmp: Vec<Vec<BackTraceNode>> = vec![vec![ BackTraceNode { times: 1, node: result } ]];
    let mut prev_idx = 0;

    loop {
        unsafe { assert_unchecked(tmp.len() >= 1); }
        let mut new_trace_node = Vec::<BackTraceNode>::new();
        let prev_ref = &tmp[prev_idx];

        for k in prev_ref {
            let compose_result = &k.node;
            let count = &k.times;

            let all_next_node = find_prev_compose_node(compose_result.clone(), selected_compose);
            for j in all_next_node {
                new_trace_node.push(BackTraceNode { times: count * j.1.1, node: j.1.0.clone() });
                for k in j.1.0.prev_gene_group {
                    if k.1 == 0 {
                        leaf.push((k.clone(), count * j.1.1));
                    }
                }
            }
        }

        let l = new_trace_node.len();

        if l != 0 {
            tmp.push(new_trace_node);
            prev_idx += 1;
        }

        if l == 0 {
            break;
        }
    }

    let mut step_count = 0usize;
    for i in (0..tmp.len()).rev() {
        if tmp[i].len() == 1 {
            output += &format!("#{}  合成{}:  概率：1/{}\n", step_count + 1, &tmp[i][0].node.gene_group, &tmp[i][0].node.probability.1);
            for genes in &tmp[i][0].node.prev_gene_group {
                output += &format!("    {}\n", genes);
            }
        } else {
            for j in 0..tmp[i].len() {
                output += &format!("#{}-{}  合成{}:  概率: 1/{}\n", step_count + 1, j + 1, &tmp[i][j].node.gene_group, &tmp[i][j].node.probability.1);
                for genes in &tmp[i][j].node.prev_gene_group {
                    output += &format!("    {}\n", genes);
                }
            }
        }

        println!();
        step_count += 1;

    }

    output
}

#[inline]
fn add_gene(genes: &mut NoHashSet<GeneGroup>) {
    println!("输入添加基因：(例XYGXYM)");
    let mut input = String::new();
    stdin().read_line(&mut input).expect("输入设备错误");
    input = input.trim().to_string();


    let splits = input.split_ascii_whitespace();
    for token in splits {  
        if token.len() != 6 {
            println!("❌ 输入有误");
            return;
        }
        let input = token.chars();
        let mut idx = 0;
        let mut g1 = [Gene::PLACEHOLDER; 6];
        for chr in input {
            if chr == 'X' || chr == 'x' {
                g1[idx] = Gene::X;
            }
            if chr == 'Y' || chr == 'y' {
                g1[idx] = Gene::Y;
            }
            if chr == 'G' || chr == 'g' {
                g1[idx] = Gene::G;
            }
            if chr == 'H' || chr == 'h' {
                g1[idx] = Gene::H;
            }
            if chr == 'W' || chr == 'w' {
                g1[idx] = Gene::W;
            }
            idx += 1;
        }

        let g = GeneGroup(g1, 0);
        if genes.contains_key(&g) {
            println!("❌ 基因已存在");
            return;
        }
        genes.insert(g, ());
    }
}

#[inline]
fn remove_gene(genes: &mut NoHashSet<GeneGroup>) {
    println!("输入删除基因：(例XYGXYM)");
    let mut input = String::new();
    stdin().read_line(&mut input).expect("❌ 输入设备错误");
    input = input.trim().to_string();
    let splits = input.split_ascii_whitespace();
    for token in splits {  
        if token.len() != 6 {
            println!("❌ 输入有误");
            return;
        }
        let input = token.chars();
        let mut idx = 0;
        let mut g1 = [Gene::PLACEHOLDER; 6];
        for chr in input {
            if chr == 'X' || chr == 'x' {
                g1[idx] = Gene::X;
            }
            if chr == 'Y' || chr == 'y' {
                g1[idx] = Gene::Y;
            }
            if chr == 'G' || chr == 'g' {
                g1[idx] = Gene::G;
            }
            if chr == 'H' || chr == 'h' {
                g1[idx] = Gene::H;
            }
            if chr == 'W' || chr == 'w' {
                g1[idx] = Gene::W;
            }
            idx += 1;
        }

        let g = GeneGroup(g1, 0);
        genes.remove(&g);
    }

}

#[inline]
fn change_target_gene() {
    println!("输入目标基因：(例3G3Y、2G1H)");
    let mut input = String::new();
    stdin().read_line(&mut input).expect("输入设备错误");
    input = input.trim().to_string();

    if input.len() % 2 != 0 || input.len() == 0 {
        println!("❌ 输入有误");
        return;
    }

    let mut input = input.chars();
    let mut idx = 0;
    
    for i in 0..5 {
        unsafe { AMOUNT_TARGET [i] = 0 };
    } 

    loop {
        let number = input.next();
        match number {
            Some(n) => {
                let chr = input.next().unwrap();
                let number = n.to_string().parse::<u8>();
                match number {
                    Ok(number) => {
                        unsafe {
                            if chr == 'G' || chr == 'g' {
                                AMOUNT_TARGET[0] = number;
                                continue;
                            }
                            if chr == 'H' || chr == 'h' {
                                AMOUNT_TARGET[1] = number;
                                continue;
                            }
                            if chr == 'Y' || chr == 'y' {
                                AMOUNT_TARGET[2] = number;
                                continue;
                            }
                            if chr == 'W' || chr == 'w' {
                                AMOUNT_TARGET[3] = number;
                                continue;
                            }
                            if chr == 'X' || chr == 'x' {
                                AMOUNT_TARGET[4] = number;
                                continue;
                            }
                            println!("❌ 输入有误");
                            return;
                        }
                    }
                    Err(_) => {
                        println!("❌ 输入有误");
                        return;
                    }
                }
            }
            None => {
                return;
            }
        }
    }
}

#[inline]
fn change_spread_limit() {
    println!("输入限制数值：");
    let mut n = "".to_string();
    stdin().read_line(&mut n).expect("❌ 输入设备错误");
    let n = n.trim();
    match n.parse::<usize>() {
        Ok(n) => {
            unsafe { SPEARD_LIMIT = n };
        }
        Err(_) => {
            println!("❌ 格式错误");
            return;
        }
    }
}

#[inline]
fn change_probability(probability: BigDecimal) -> Option<BigDecimal> {
    print!("输入概率：1/");
    stdout().flush();
    let mut n = "".to_string();
    stdin().read_line(&mut n).expect("❌ 输入设备错误");
    let n = n.trim();
    match n.parse::<usize>() {
        Ok(n) => {
            return Some(BigDecimal::from_usize(n).unwrap());
        }
        Err(_) => {
            println!("❌ 格式错误");
            return None;
        }
    }
}

fn main() {
    // let mut genes_vec: Vec<GeneGroup> = genes!{ "GYYHYY", "GYYYGY", "GYYYYY", "GGGHYX", "XYGHYW", "GYXYXY", "XYHGGY", };
    let mut genes_vec = nohashset!(50);
    let mut probability_filter = BigDecimal::from_usize(1000).unwrap();

    loop {
        let mut global_lowest_loss: Option<ComposeResult> = None;
        let mut already_compose_collection = NoHashSet::<ComposeResult>::with_capacity_and_hasher(200, BuildNoHashHasher::default());
        let mut selected_composed = NoHashSet::<ComposeResult>::with_capacity_and_hasher(100, BuildNoHashHasher::default());
        println!("\n\n\n\n\n\n\n\n\n\n\n\n🎯 目标基因：{}", display_target());
        print!("🧬 当前基因：");
        for g in &genes_vec {
            print!("{} ", g.0.display())
        }
        println!("\n🚫 概率过滤器: >=1/{}", probability_filter );
        println!("🚫 最长搜索链: <{}", unsafe { SPEARD_LIMIT });
        println!();
        println!("1. 🔵 添加基因");
        println!("2. 🔵 删除基因");
        println!("3. 🔵 修改目标基因");
        println!("4. 🔵 清除所有基因");
        println!("5. 🟡 修改搜索链限制");
        println!("6. 🟡 修改概率过滤器");
        println!("9. 开始");
        println!("0. 退出");

        print!("\n> ");
        stdout().flush();

        let mut code = String::new();
        io::stdin().read_line(&mut code).expect("输入出错");
        code = code.trim().to_string();

        if code == "0" {
            exit(0);
        }

        if code == "1" {
            add_gene(&mut genes_vec);
        }

        if code == "2" {
            remove_gene(&mut genes_vec);
        }

        if code == "3" {
            change_target_gene();
        }

        if code == "4" {
            genes_vec = nohashset!(50);
        }

        if code == "5" {
            change_spread_limit();
        }

        if code == "6" {
            let r = change_probability(probability_filter);
            probability_filter = match r {
                Some(r) => r,
                None => BigDecimal::from_usize(1000).unwrap(),
            }
        }

        if code == "9" {
            if genes_vec.len() == 0 {
                println!("❌ 基因组不足");
            } else {
        
                let mut genes_vec = genes_vec.clone();

                let start_time = SystemTime::now();
                let (final_compose, partial_lowest) = 
                    predict(&mut genes_vec, &mut already_compose_collection, &mut selected_composed, &probability_filter);
                
                match final_compose {
                    Some(a) => {
                        println!("🟢 ✔ 找到路径  🟢\n{}", display_backtrace_path(a, &mut selected_composed));
                    }
                    None => {
                        println!("🔴 ❌ 无路径  🔴");

                        match partial_lowest {
                            Some(c) => { 
                                println!("但是找到一个推荐路径：\n{}", display_backtrace_path(c, &mut selected_composed));
                            }
                            None => { }
                        }
                    }
                }
                let end_time = start_time.elapsed();
                match end_time {
                    Ok(t) => {
                        println!("用时: {}ms", t.as_millis());
                    }
                    Err(_) => { },
                }
            }
        }

        println!("\n按下回车继续.");
        let mut g = "".to_string();
        stdin().read_line(&mut g).expect("");
        unsafe { GLOBAL_SIGN = 1 };
    }
}