(require (lib "defmacro.ss"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Top level functions for regression or numerical evaluation of belief

;; get the formula that is the regression of (BEL wff) after acts
(define-macro (regr-bel wff acts . args) `(regress '(BEL ,wff) ',acts))

;; get an estimated degree of belief in wff after acts with some timing info
(define-macro (eval-bel wff acts . args) 
  (set! error-sigma (if (null? args) default-error-sigma (car args)))
  `(begin 
     (printf "Regression ")
     (let ((expr (time (regr-bel ,wff ,acts))))
       (printf "Estimation of belief with sigma=~a\n" ,error-sigma)
       (time (eval (make-monte expr))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Handle the BAT spec.  See README for format.  
;;;  Global vars: all-fluents, ini-p-expr, l-exprs, act-alts, ss-exprs, 

(define-macro (define-fluents . flus) `(define all-fluents ',flus))
(define-macro (define-ini-p-expr expr) `(define ini-p-expr ,expr))
(define-macro (define-l-exprs . pairs) `(define l-exprs (hasheq ,@(pp pairs))))
(define-macro (define-alts . pairs) `(define act-alts (hasheq ,@(pp pairs))))
(define ss-exprs (make-hasheq))
(define-macro (define-ss-exprs flu . pairs)
  `(hash-set! ss-exprs ',flu (hasheq ,@(pp pairs))))

;; convert list ... (act . args) expr ... to ... 'act (lambda args expr) ...
(define-for-syntax (pp pairs)
  (if (null? pairs) '()
      (cons `',(caar pairs) 
            (cons `(lambda ,(cdar pairs) ,(cadr pairs)) (pp (cddr pairs))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Regression 

;; T[t,acts] and R[phi,acts]. Note: the interpretation of the action sequence. E.g., (sonar 5)(nfwd 2 2) means sonar was done first, so we reverse this order and regression begins with the last action performed
(define (regress e revacts)
  (define acts (reverse revacts))
  (cond 
   [(memq e all-fluents)   ; e is a fluent
    (if (null? acts) e (regress (apply-ss-axiom e (car acts)) (cdr acts)))]
   [(not (list? e)) e]     ; e is a number or constant
   [(eq? (car e) 'BEL)     ; e is (bel wff)
    (list '/ (make-integral (cadr e) acts) (make-integral #t acts))]
   [else (cons (car e) (map (lambda (arg) (regress arg acts)) (cdr e)))]))

;; expression for successor value of fluent after act
(define (apply-ss-axiom flu act)
  (let ((ss-hash (hash-ref ss-exprs flu #f)))
    (if ss-hash
        (let ((ss-expr (hash-ref ss-hash (car act) #f)))
          (if ss-expr (apply ss-expr (cdr act)) flu))
        flu)))

;; expression for likelihood of act
(define (apply-l-axiom act)
  (let ((l-expr (hash-ref l-exprs (car act) #f)))
    (if l-expr (apply l-expr (cdr act)) 1.0)))

;; the integral/summation of the density function for states where wff is true
(define (make-integral wff acts)
  (let-values (([alts vars] (gen-alts (reverse acts) '() '())))
    `(INTEGRATE ,(append all-fluents (reverse vars))
                    ,(product-expr ini-p-expr (regress-dense wff alts)))))

;; given acts, return alts and a list of new prob variables required
(define (gen-alts acts vars alts)
  (if (null? acts) (values alts vars)
      (let ((altfn (hash-ref act-alts (caar acts) #f)))
        (if (not altfn) (gen-alts (cdr acts) vars (cons (car acts) alts))
            (let* ((var (gensym 'z)) (alt ((apply altfn (cdar acts)) var)))
              (gen-alts (cdr acts) (cons var vars) (cons alt alts)))))))

;; T[P(wff),acts]
(define (regress-dense wff acts)
  (if (null? acts) (if (eq? wff #t) 1.0 `(if ,wff 1.0 0.0))
      (let ((lval (regress (apply-l-axiom (car acts)) (cdr acts)))
            (dval (regress-dense (regress wff (list (car acts))) (cdr acts))))
        (product-expr lval dval))))

;; an expression that is equal to the product of exprs, but in unnested form
(define (product-expr . exprs)
  (let loop ((exprs exprs))
    (if (null? exprs) 1.0
        (let ((e1 (car exprs)) (e2 (loop (cdr exprs))))
          (cond 
           [(equal? e1 1.0) e2]    [(equal? e2 1.0) e1] 
           [(and (list? e1) (eq? (car e1) '*))
            (if (and (list? e2) (eq? (car e2) '*)) `(* ,@(cdr e1) ,@(cdr e2))
                `(* ,@(cdr e1) ,e2))]
           [(and (list? e2) (eq? (car e2) '*)) `(* ,e1 ,@(cdr e2))]
           [else `(* ,e1 ,e2)])))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   Known probability distributions (for convenience)

;; normal distribution of x with mean mu and variance var
(define (GAUSSIAN x mu var)
  (/ (exp (- (/ (* (- x mu) (- x mu)) (* 2.0 var)))) (sqrt (* 2.0 var pi))))

;; uniform distribution of x with lo and high
(define (UNIFORM x low high)
  (if (and (> x low) (< x high)) (/ 1.0 (- high low)) 0.0))

;; discrete distribution with values vi from v1 p1 ... vn pn, where sum of pi=1
(define (DISCRETE x . args)
  (let loop ((args args))
    (if (null? args) 0.0
        (if (eq? (car args) x) (cadr args) (loop (cddr args))))))


(define (BINARY x p) (DISCRETE x #t p #f (- 1 p)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Convert expr with (INTEGRATE-SUM ...) into the (MONTE-CARLO ...) form
;;;  for numerical integration below, or produce an error otherwise

(define (make-monte e)
  (cond
   [(not (list? e)) e]
   [(eq? (car e) 'INTEGRATE)
    (let ((vars (cadr e)) (expr (caddr e)))
      (if (and (list? expr) (eq? (car expr) '*))
          (let-values (([vspecs mults] (get-vspecs vars (cdr expr))))
            `(MONTE-CARLO ,vspecs ,(apply product-expr mults)))
          (error "error: integrand is not a product" expr)))]
   [else (cons (car e) (map make-monte (cdr e)))]))

;; return vspecs for vars from mults (exprs to be multiplied) and mults left
(define (get-vspecs vars mults)
  (let ((dists (map (lambda (v) (get-var-dist v mults)) vars)))
    (values (map (lambda (v d) (list v (dist->gen v d))) vars dists)
            (filter (lambda (m) (not (memq m dists))) mults))))

;; first expression in mults that gives a known distribution for v
(define (get-var-dist v mults)
  (ormap (lambda (m) (and (list? m) (list? (cdr m)) (eq? (cadr m) v)
                          (hash-ref dist-trans (car m) #f) m))
         mults))

;; replace distribution expression by random value generator expression
(define (dist->gen var expr)
  (if (not expr) (error "error: no distribution expression for variable" var)
      (cons (hash-ref dist-trans (car expr)) (cddr expr))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Random value generators for the three known distributions

;; mapping between distributions and generators for the distribution
(define dist-trans (hasheq 'GAUSSIAN 'GAUSSIAN-RAND 
                           'UNIFORM 'UNIFORM-RAND 
                           'DISCRETE 'DISCRETE-RAND))

;; box muller transform for generating Gaussian numbers
(define (GAUSSIAN-RAND mu sigma) 
  (let loop ()
    (let* ((z1 (- (* 2 (random)) 1)) (z2 (- (* 2 (random) 1)))
           (rsq (+ (* z1 z1) (* z2 z2))))
      (if (> rsq 1) (loop)
          (let ((x1 (* z1 (sqrt (/ (* -2 (log rsq)) rsq)))))
            (+ (* x1 sigma) mu))))))

;; linear transform for generating uniform numbers
(define (UNIFORM-RAND low high) (+ low (* (- high low) (random))))

;; generating discrete values vi for args as in DISCRETE above
(define (DISCRETE-RAND . args)
  (let ((r (random)))
    (let loop ((args args) (sum 0))
      (let ((next (+ sum (cadr args))))
        (if (<= r next) (car args) (loop (cddr args) next))))))

;; generating binary values with given prob for #t 
(define (BINARY-RAND p) (DISCRETE-RAND #t p #f (- 1 p)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Numerical integration

(define-for-syntax default-error-sigma .01)
(define-for-syntax error-sigma default-error-sigma)

;; integral [-inf,inf] of f(x)=expr*q(x) where q is prob distr given by vspecs
(define-macro (MONTE-CARLO vspecs expr . args)
  (let ((sigma (if (null? args) error-sigma (car args))))
    `(begin
       (eprintf "Integration over ~a variables: N is " ,(length vspecs))
       (monte-sample (lambda () (let* ,vspecs ,expr)) ,sigma))))

;; get average of random f values, checking every 1000 for estimated error
(define (monte-sample f sigma)
  (define variance (* sigma sigma))
  (let loop ((n 0) (sum 0) (sumsq 0))
    (if (and (> n 1) (< (- (* n sumsq) (* sum sum)) (* variance (- n 1) n n)))
        (begin (eprintf "~a x 10^6\n" (/ n 1000000.0)) (/ sum n))
        (let loop2 ((m 1000) (s2 0) (sq2 0))
          (if (= m 0) (loop (+ n 1000) (+ sum s2) (+ sumsq sq2))
              (let ((y (f))) (loop2 (- m 1) (+ s2 y) (+ sq2 (* y y)))))))))

