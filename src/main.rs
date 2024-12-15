use std::{borrow::Borrow, cell::RefCell, collections::{HashMap, HashSet}, fmt::Display, fs::exists, hint::assert_unchecked, io::{self, stdin, stdout, Write}, process::{exit, id}, rc::Rc, sync::LazyLock, usize, vec};

use bigdecimal::{BigDecimal, FromPrimitive};
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
            sign: 0,
            loss: RefCell::new(u8::MAX),
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
fn select2compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>, probability: &BigDecimal) -> Vec<ComposeResult> {
    // already_collection only for check repeat item
    let mut output = Vec::new();
    let l = exist_gene.len();
    for i in 0..l {
        for j in i..l {
            if i != j {
                let v = cl! { exist_gene[i], exist_gene[j], };
                let results = compose(v);
                let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                output.append(&mut added);
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select3compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>, probability: &BigDecimal) -> Vec<ComposeResult> {
    // already_collection only for check repeat item
    let mut output = Vec::new();
    let l = exist_gene.len();
    for i in 0..l {
        for j in i..l {
            for k in j..l {
                if !(i == j && i == k) {
                    let v = cl! { exist_gene[i], exist_gene[j], exist_gene[k], };
                    let results = compose(v);
                    let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                    output.append(&mut added);
                }
            }
        }
    }
    output
}

/// Output is the added items for already_collection
fn select4compose(exist_gene: &Vec<GeneGroup>, already_collection: &Vec<ComposeResult>, probability: &BigDecimal) -> Vec<ComposeResult> {
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
                        let mut added = vec_result_non_exists(&results, &output, &already_collection, probability);
                        output.append(&mut added);
                    }
                }
            }
        }
    }
    output
}
fn vec_result_non_exists(result: &Vec<ComposeResult>, v0: &Vec<ComposeResult>, v1: &Vec<ComposeResult>, probability: &BigDecimal) -> Vec<ComposeResult> {
    let mut output = Vec::new();
    for r in result {
        let probability_denominator = &r.probability.1;
        if !result_exists(r, v0) && !result_exists(r, v1) && (probability_denominator <= probability) {
            output.push(r.clone());
        }
    }
    output
}

fn vec_result_non_selected(result: Vec<ComposeResult>, selected_collection: &Vec<ComposeResult>, probability: &BigDecimal) -> Vec<ComposeResult> {
    let mut output = Vec::new();
    if selected_collection.len() == 0 {
        return result;
    }
    'outter: for i in result {
        for selected in selected_collection {
            if i != *selected && (&i.probability.1 <= probability) {
                output.push(i);
                continue 'outter;
            }
        }
    }
    output
}

#[inline]
fn result_exists(result: &ComposeResult, v: &Vec<ComposeResult>) -> bool {
    let mut exists = false;
    for n in v {
        if n == result {
            exists = true;
        }
    }
    exists
}

fn one_step_compose_predict<'a>(exists_gene: &Vec<GeneGroup>, already_compose_collection: &Vec<ComposeResult>, selected_collection: &Vec<ComposeResult>, probability: &BigDecimal) 
    -> Vec<ComposeResult> {
    let mut a1 = select2compose(&exists_gene, &already_compose_collection, probability);
    let mut a2 = select3compose(&exists_gene, &already_compose_collection, probability);
    let mut a3 = select4compose(&exists_gene, &already_compose_collection, probability);


    let mut a2 = vec_result_non_exists(&a2, &a1, &already_compose_collection, probability);
    a1.append(&mut a2);
    let mut a3 = vec_result_non_exists(&a3, &a1, &already_compose_collection, probability);
    a1.append(&mut a3);
    a1 = vec_result_non_selected(a1, selected_collection, probability);

    for j in &mut a1 {
        // SAFTY: Use static data with single thread is safe. No any data conflition.
        j.sign = unsafe { GLOBAL_SIGN };
        j.gene_group.1 = unsafe { GLOBAL_SIGN };
        unsafe { GLOBAL_SIGN += 1 };
    }

    evaluate_loss(&mut a1);
    a1
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

fn predict(exist_gene: &mut Vec<GeneGroup>, composed_collection: &mut Vec<ComposeResult>, selected_collection: &mut Vec<ComposeResult>, probability: &BigDecimal) -> (Option<ComposeResult>, Option<ComposeResult>) {
    let mut generation = 1usize;
    let mut partial_lowest: Option<ComposeResult> = None;
    loop {
        let mut added_compose: Vec<ComposeResult> = one_step_compose_predict(&exist_gene, &composed_collection, &selected_collection, probability);
        composed_collection.append(&mut added_compose); 
        if composed_collection.len() == 0 || exist_gene.len() > unsafe { SPEARD_LIMIT } {
            return (None, partial_lowest);
        }

        unsafe { assert_unchecked((&composed_collection).len() > 0); }
        // Find lowest loss
        let mut lowest_loss: &ComposeResult = &composed_collection[0];
        let mut lowest_idx = 0usize;
        for k in 0..composed_collection.len() {
            if *(&composed_collection[k]).loss.borrow() < *lowest_loss.loss.borrow() {
                lowest_loss = &composed_collection[k];
                lowest_idx = k;
            }
        }

        match partial_lowest {
            Some(ref c) => {
                if lowest_loss.loss < c.loss {
                    partial_lowest = Some(lowest_loss.clone());
                }
            }
            None => {
                partial_lowest = Some(lowest_loss.clone());
            }
        }

        if *lowest_loss.loss.borrow() != 0 {
            exist_gene.push(lowest_loss.gene_group.clone());
            selected_collection.push(lowest_loss.clone());
            composed_collection.remove(lowest_idx);
        } else {
            return (Some(lowest_loss.clone()), partial_lowest);
        }
    }
    generation += 1;
}

macro_rules! genes {
    ($($s: expr,)+) => {
        vec![$(GeneGroup::from_str($s).unwrap(),)+]
    };
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

fn find_prev_compose_node(node: ComposeResult, selected_collection: &Vec<ComposeResult>) -> HashMap<usize, (ComposeResult, usize)> {
    let mut v = HashMap::<usize, (ComposeResult, usize)>::new();
    for k in node.prev_gene_group {
        if k.1 != 0 {
            if !v.contains_key(&k.1) {
                // Find compose result
                for c in selected_collection {
                    if c.sign == k.1 {
                        v.insert(k.1, (c.clone(), 1));
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

fn display_backtrace_path(result: ComposeResult, selected_compose: &Vec<ComposeResult>) -> String {

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
fn add_gene(genes: &mut Vec<GeneGroup>) {
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
        for i in 0..genes.len() {
            if g == genes[i] {
                println!("❌ 基因已存在");
                return;
            }
        }
        genes.push(g);
    }

}

#[inline]
fn remove_gene(genes: &mut Vec<GeneGroup>) {
    println!("输入删除基因：(例XYGXYM)");
    let mut input = String::new();
    stdin().read_line(&mut input).expect("❌ 输入设备错误");
    input = input.trim().to_string();
    let splits = input.split_ascii_whitespace();
    'outter: for token in splits {  
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
        for i in 0..genes.len() {
            if genes[i] == g {
                genes.remove(i);
                continue 'outter;
            }
        }
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
    let mut genes_vec = Vec::<GeneGroup>::new();
    let mut probability_filter = BigDecimal::from_usize(1000).unwrap();

    loop {
        let mut global_lowest_loss: Option<ComposeResult> = None;
        let mut already_compose_collection = Vec::new();
        let mut selected_composed = Vec::new();
        println!("\n\n\n\n\n\n\n\n\n\n\n\n🎯 目标基因：{}", display_target());
        print!("🧬 当前基因：");
        for g in &genes_vec {
            print!("{} ", g.display())
        }
        println!("\n🚫 概率过滤器: >=1/{}", probability_filter );
        println!("🚫 最长搜索链: <{}", unsafe { SPEARD_LIMIT });
        println!();
        println!("1. 🟢 添加基因");
        println!("2. 🟢 删除基因");
        println!("3. 🟢 修改目标基因");
        println!("4. 🟢 清除所有基因");
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
            change_spread_limit();
        }

        if code == "5" {
            let r = change_probability(probability_filter);
            probability_filter = match r {
                Some(r) => r,
                None => BigDecimal::from_usize(1000).unwrap(),
            }
        }

        if code == "9" {
            genes_vec = vec![];
        }

        if code == "6" {
            if genes_vec.len() == 0 {
                println!("❌ 基因组不足");
            } else {
        
            let mut genes_vec = genes_vec.clone();
            let (final_compose, partial_lowest) = predict(&mut genes_vec, &mut already_compose_collection, &mut selected_composed, &probability_filter);
            
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
            }
        }

        println!("\n按下回车继续.");
        let mut g = "".to_string();
        stdin().read_line(&mut g).expect("");
        unsafe { GLOBAL_SIGN = 1 };
    }
}