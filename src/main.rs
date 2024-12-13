use std::{collections::HashSet, hint::assert_unchecked, vec};

use fraction::Fraction;

mod fraction;

macro_rules! cl {
    ($($name: expr,)+) => {
        vec! [$($name.clone(),)+]
    };
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

#[derive(Debug)]
#[derive(Clone)]
struct ComposeResult {
    sign: usize,
    gene_group: GeneGroup,
    probability: Fraction,
    prev_gene_group: Vec<GeneGroup>,
}

impl PartialEq for ComposeResult {
    fn eq(&self, other: &Self) -> bool {
        self.gene_group == other.gene_group
    }
}

impl PartialEq for GeneGroup {
    fn eq(&self, other: &Self) -> bool {
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
}

impl Clone for GeneGroup {
    fn clone(&self) -> Self {
        Self(self.0, self.1)
    }
}

static mut GLOBAL_SIGN: usize = 1;

fn compose(plants: Vec<GeneGroup>) -> Vec<ComposeResult> {
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
    let mut composes = Vec::<ComposeResult>::new();
    for k in expanded {
        let c = ComposeResult {
            // SAFTY: Use static data with single thread is safe. No any data conflition.
            sign: 0,
            gene_group: k.0,
            probability: k.1,
            prev_gene_group: plants.clone(),
        };
        composes.push(c);
    }

    composes
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

/// Output is the added items for already_collection
fn select2compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>) -> Vec<ComposeResult> {
    // already_collection only for check repeat item
    let mut output = Vec::new();
    let l = exist_gene.len();
    for i in 0..l {
        for j in i..l {
            if i != j {
                let v = cl! { exist_gene[i], exist_gene[j], };
                let results = compose(v);
                let mut added = vec_result_non_exists(&results, &output, &already_collection);
                output.append(&mut added);
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select3compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>) -> Vec<ComposeResult> {
    // already_collection only for check repeat item
    let mut output = Vec::new();
    let l = exist_gene.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                if !(i == j && i == k) {
                    let v = cl! { exist_gene[i], exist_gene[j], exist_gene[k], };
                    let results = compose(v);
                    let mut added = vec_result_non_exists(&results, &output, &already_collection);
                    output.append(&mut added);
                }
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select4compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>) -> Vec<ComposeResult> {
    // already_collection only for check repeat item
    let mut output = Vec::new();
    let l = exist_gene.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                for z in k..l {
                    if !(i == j && i == k && i == z) {
                        let v = cl! { exist_gene[i], exist_gene[j], exist_gene[k], };
                        let results = compose(v);
                        let mut added = vec_result_non_exists(&results, &output, &already_collection);
                        output.append(&mut added);
                    }
                }
            }
        }
    }
    output
}
fn vec_result_non_exists(result: &Vec<ComposeResult>, v0: &Vec<ComposeResult>, v1: &Vec<ComposeResult>) -> Vec<ComposeResult> {
    let mut output = Vec::new();
    for r in result {
        if !result_exists(r, v0) && !result_exists(r, v1) {
            output.push(r.clone());
        }
    }
    output
}

fn result_exists(result: &ComposeResult, v: &Vec<ComposeResult>) -> bool {
    let mut exists = false;
    for n in v {
        if n == result {
            exists = true;
        }
    }
    exists
}

fn one_step_compose_predict(exists_gene: &Vec<GeneGroup>, already_compose_collection: &Vec<ComposeResult>) -> Vec<ComposeResult> {
    let mut a1 = select2compose(&exists_gene, &already_compose_collection);
    let mut a2 = select3compose(&exists_gene, &already_compose_collection);
    let mut a3 = select4compose(&exists_gene, &already_compose_collection);

    let mut a2 = vec_result_non_exists(&a2, &a1, &already_compose_collection);
    a1.append(&mut a2);
    let mut a3 = vec_result_non_exists(&a3, &a1, &already_compose_collection);
    a1.append(&mut a3);
    for j in &mut a1 {
        // SAFTY: Use static data with single thread is safe. No any data conflition.
        j.sign = unsafe { GLOBAL_SIGN };
        j.gene_group.1 = unsafe { GLOBAL_SIGN };
        unsafe { GLOBAL_SIGN += 1 };
    }
    a1
}

macro_rules! genes {
    ($($s: expr,)+) => {
        vec![$(GeneGroup::from_str($s).unwrap(),)+]
    };
}


fn main() {
    let genes_vec = genes! { "WGYXYW", "XYYGYG", "XGWYGY", "XYWWGY", "GYYHWX", "XYYXYX", };
    let already_compose_collection = Vec::new();

    let output = one_step_compose_predict(&genes_vec, &already_compose_collection);

    for j in &output {
        println!("{:?} -> {:?} #{:?}", &j.gene_group, &j.probability, &j.sign);
    }

}