(include "../regress.scm")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  bat only with horizontal distance 

(define-fluents h)

;; initially, all four fluents are independent variables
(define-ini-p-expr 
  `(UNIFORM h 2 12)           ; h is uniform in range 2 12 
   )

;; for each noisy action, a function that returns a variant expression
(define-alts
   (nfwd x y) (lambda (z) `(nfwd ,x ,z)))

;; for each action, an expression for its likelihood (default = 1)
(define-l-exprs
   (sonar z) `(GAUSSIAN ,z h 4.0)            ; true value + noise
   (nfwd x y) `(GAUSSIAN ,y ,x 1.0))   ; noisy arg slightly favors #t

;; for each action, an expression for the next value of h (default = as is)
(define-ss-exprs h
   (fwd z) `(max 0 (- h ,z))
   (nfwd x z) `(max 0 (- h ,z)) )

;; for each action, an expression for the next value of light (default = as is)

