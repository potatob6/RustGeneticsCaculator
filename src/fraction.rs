use std::{fmt::Debug, hint::assert_unchecked, ops::{Add, Mul}};

use bigdecimal::{BigDecimal, FromPrimitive, Zero};

#[derive(Clone)]
pub struct Fraction(BigDecimal, BigDecimal);

impl Default for Fraction {
    fn default() -> Self {
        Self(BigDecimal::zero(), BigDecimal::zero())
    }
}

impl Fraction {
    #[inline(always)]
    pub fn reduce(&self) -> Self {
        let gcd = gcd(&self.0, &self.1);
        Self(&self.0 / &gcd, &self.1 / &gcd)
    }

    #[inline(always)]
    pub fn from_u32(numerator: u32, denominator: u32) -> Self {
        Self(BigDecimal::from_u32(numerator).unwrap(), BigDecimal::from_u32(denominator).unwrap())
    }
}

impl Debug for Fraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = format!("{}/{}", self.0, self.1);
        f.write_str(&s)
    }
}

fn gcd(x: &BigDecimal, y: &BigDecimal) -> BigDecimal {
    unsafe { assert_unchecked(x >= &BigDecimal::zero()); }
    unsafe { assert_unchecked(y >= &BigDecimal::zero()); }
    if x.eq(y) {
        return x.clone();
    }

    let mut r1 = if x >= y { x } else { y };
    let mut r2 = if x >= y { y } else { x };

    if r2.eq(&BigDecimal::zero()) {
        return BigDecimal::zero();
    }

    let mut v_stack = Vec::<(BigDecimal, BigDecimal)>::new();

    while !r2.eq(&BigDecimal::zero()) {
        v_stack.push((r2.clone(), r1 % r2));
        let last = v_stack.last();
        unsafe { assert_unchecked(last.is_some()); }
        let last = last.unwrap();
        r1 = &last.0;
        r2 = &last.1;
    }

    r1.clone()

}


impl Add<&Fraction> for Fraction {
    #[inline(always)]
    fn add(self, rhs: &Self) -> Self::Output {
        &self + rhs
    }
    
    type Output = Self;
}

impl<'a, 'b> Add<&'a Fraction> for &'b Fraction {
    type Output = Fraction;

    #[inline(always)]
    fn add(self, rhs: &'a Fraction) -> Self::Output {
        let dem_gcd: BigDecimal = gcd(&self.1, &rhs.1);
        let lcm = &dem_gcd * (&self.1 / &dem_gcd) * (&rhs.1 / &dem_gcd);
        let f = Fraction(&self.0 * &lcm / &self.1 + &rhs.0 * &lcm / &rhs.1, lcm);
        f.reduce()
    }
}

impl Mul<&Fraction> for Fraction {
    #[inline(always)]
    fn mul(self, rhs: &Self) -> Self::Output {
        let f = &self * rhs;
        f.reduce()
    }
    
    type Output = Self;
}

impl<'a, 'b> Mul<&'a Fraction> for &'b Fraction {
    type Output = Fraction;

    #[inline(always)]
    fn mul(self, rhs: &'a Fraction) -> Self::Output {
        let k = Fraction(&self.0 * &rhs.0, &self.1 * &rhs.1);
        k.reduce()
    }
}

