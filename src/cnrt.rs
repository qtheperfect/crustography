/**
 * Infor:  The Chinese Remainder Theorem Instructional Program.
 * Author: Linghuei Guo
 * Contact: PGP - 6819D81B0971C2C4
 * Permission:   All Rights Reserved
 */
pub mod cnrt {

    use num_bigint::BigInt;

    fn into_mod(x: &BigInt, m_r: &BigInt) -> BigInt {
	(x % m_r + m_r) % m_r
    }

    fn bezout_marked(x: &BigInt, y: &BigInt, mark: i32) -> (BigInt, BigInt, BigInt) {
	/** Find (u, v, c) s.t.
	 *   u > 0, v <= 0, coprime(u, v),
	 *   ux + vy = c */
	if *y == BigInt::ZERO {
	    (BigInt::from(1), BigInt::from(mark), x.clone())
	}
	else {
	    let x1 = x % y;
	    let k = x / y;
	    let (mut v, u, c) = bezout_marked(y, &x1, 1 - mark);
	    v = &v - &u * &k;
	    (u, v, c)
	}
    } // end fn

    pub fn bezout_recursive(x: &BigInt, y: &BigInt) -> (BigInt, BigInt, BigInt) {
	bezout_marked(x, y, 0)
    }

    pub fn bezout(x: &BigInt, y: &BigInt) -> (BigInt, BigInt, BigInt) {
	fn looper(mut x1: BigInt, x2: BigInt, mut s11: BigInt, mut s21: BigInt, mut s12: BigInt, mut s22: BigInt)
		  -> (BigInt, BigInt, BigInt) {
	    if x2 == BigInt::ZERO {
		(s11, s21, x1)
	    } else {
		let k = &x1 / &x2;
		x1 = &x1 % &x2;
		s11 = &s11 - &k * &s12;
		s21 = &s21 - &k * &s22;
		looper (x2, x1, s12, s22, s11, s21)
	    }
	}
	looper(x.clone(), y.clone(), BigInt::from(1), BigInt::from(0), BigInt::from(0), BigInt::from(1))
    }

    pub fn test_bezout() {
	let x = "828342".parse::<BigInt>().unwrap();
	let y = "423512344114231524".parse::<BigInt>().unwrap();
	let (u, v, c) = bezout(&x, &y);
	println!("{} * {} + {} * {} = {} == (gcd: {})",
		 &u, &x, &v, &y,
		 &u * &x + &v * &y, &c
	)
    }

    #[derive(Clone)]
    pub struct RemainderValue {
	r: BigInt,
	m: BigInt
    }

    /**
     * Find the Chinese Remainder Problem Solution r so that
     *   r == r1 (mod m1) && r == r2 (mod m2) where:
     *    gcd(m1, m2) == 1 */
    fn find_cnrt(m1: &BigInt, m2: &BigInt, r1: &BigInt, r2: &BigInt) -> RemainderValue {
	let (u, v, c) = bezout(m1, m2);
	//let ans = r1 + (r2 - r1) * &u * m1;
	let ans = r2 * &u * m1 + r1 * &v * m2;
	let mod12 = m1 * m2;
	RemainderValue {
	    r: into_mod(&ans, &mod12), 
	    m: mod12
	}
    }


    impl RemainderValue {

	pub fn new() -> RemainderValue {
	    RemainderValue {
		r: BigInt::from(0),
		m: BigInt::from(1),
	    }
	}

	pub fn make(n: &BigInt, m: &BigInt) -> RemainderValue {
	    RemainderValue {
		r: into_mod(n, m),
		m: m.clone()
	    }
	}

	pub fn extend(self, m1: &BigInt, r1: &BigInt) -> RemainderValue {
	    // gcd(self.m, m1) == 1 presumed.
	    find_cnrt(&self.m, m1, &self.r, r1)
	}

	
	pub fn merge(self, result2: &RemainderValue) -> RemainderValue {
	    /* Merge the results of self and result2 into a remainder
	         of a larger modular: lcm(self.m, result2.m)
	         even if gcd(self.m, result2.m) != 1 */ 
	    let (_, _, c) = bezout(&self.m, &result2.m);
	    if BigInt::from(1) == c { // ??? &c == BigInt::from(1)
		self.extend(&result2.m, &result2.r)
	    }
	    else {
		let self_prevail_factors = &self.m / &c;

		// Get rid of the prime factors in result2.m who have HIGHER orders as factors of self.m:
		let obj2 = result2.clone()
		    .exclude(&self_prevail_factors); 

		// Get rid of prime factors of self.m who have SAME OR HIGHER orders as factors of result2.m:
		self.exclude(&obj2.m)
		// and then extend into the modulos lcm(self.m, result2.m)
		    .extend(&obj2.m, &obj2.r)
	    }
	}

	pub fn exclude(self, target: &BigInt) -> RemainderValue {
	    let (_, _, c) = bezout(&self.m, target);
	    if c != BigInt::from(1) {
		let m1 = &self.m / &c;
		let ans = RemainderValue {
		    r: &self.r % &m1,
		    m: m1
		};
		ans.exclude(&c)
	    }
	    else {
		self
	    }
	}

	pub fn verify(&self, n: &BigInt) -> bool {
	    println!(" {{remainder: {}}} == {}( <{}> mod {} )",
		     &self.r,
		     n % &self.m,
		     n,
		     &self.m
	    );
	    n % &self.m == self.r
	}
    }
    
    pub fn test() {
	let num = "19122025".parse::<BigInt>().unwrap();
	let testers = vec![32, 12, 28, 77, 93, 121, 17, 711];
	println!("\n\tReconstructing the number {} from its remainders on these integers:\n\t {:?}", num, testers);

	let result = testers.iter()
	    .map(|&t| {BigInt::from(t)})
	    .map(|m| {RemainderValue::make(&num, &m)}) 
	    .fold(
		RemainderValue::new(),
		|r, cr1| {
		    print!("\n => merging:  ");
		    cr1.verify(&num);
		    let result = r.merge(&cr1);
		    result.verify(&num);
		    result
		});
    }
}

