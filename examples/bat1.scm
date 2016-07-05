(include "../regress.scm")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; This is an example BAT with four fluents:
;;; h is distance to the wall (continuous)
;;; v is vertical position (continuous)
;;; light is on or off (discrete)
;;; temp is cold, cool, tepid, warm, or hot (discrete)
;;; The physical actions are
;;; (fwd x) is reduce h by x to a minimum of 0.0
;;; (nfwd x y) is a noisy version of (fwd x): 
;;; y is real prob variable that determines amount moved given x
;;; (toggle) is toggle the light setting
;;; (ntoggle b) is a noisy version of (toggle):
;;; b is boolean prob variable that determines if toggle succeeds
;;; The sensing actions are
;;; (sonar z) read z as h state
;;; (sense-light z) read z as the light state
;;; (sense-temp z) read z as the temp state

(define-fluents h v light temp)

;; initially, all four fluents are independent variables
(define-ini-p-expr 
`(* (UNIFORM h 2 12) ; h is uniform in range 2 12 
(GAUSSIAN v 0 1) ; v is normal with mean 0 dev 1
(DISCRETE light 'on 1.0) ; light is known to be on
(DISCRETE temp 'cold .1 'cool .2 'tepid .4 'warm .2 'hot .1) ))

;; for each noisy action, a function that returns a variant expression
(define-alts
(nfwd x y) (lambda (z) `(nfwd ,x ,z)) 
(ntoggle y) (lambda (z) `(ntoggle ,z)))

;; for each action, an expression for its likelihood (default = 1)
(define-l-exprs
(sonar z) `(GAUSSIAN ,z h 4.0) ; true value + noise (sense-light z)(DISCRETE light ,z 1.0) ; exact sensing
(sense-temp z) `(if (eq? temp ,z) .6 (if (one-away? temp ,z) .2 0.0)) 
(nfwd x y) `(GAUSSIAN ,y ,x 1.0) ; noisy arg is normal around x
(ntoggle z) `(DISCRETE ,z #t .6 #f .4)) ; noisy arg slightly favors #t

;; for each action, an expression for the next value of h (default = as is)
(define-ss-exprs h
(fwd z) `(max 0 (- h ,z)) 
(nfwd x z) `(max 0 (- h ,z))  )

;; for each action, an expression for the next value of light (default = as is)
(define-ss-exprs light
(toggle) `(opposite light) 
(ntoggle z) `(if ,z (opposite light) light))

;; utility function: the opposite value of a light
(define (opposite x) (if (eq? x 'on) 'off 'on))

;; utility function: are t1 and t2 one step away in temp scale?
(define (one-away? t1 t2) 
(case t1
[(cold) (memq t2 '(cool))]
[(cool) (memq t2 '(cold tepid))]
[(tepid) (memq t2 '(cool warm))]
[(warm) (memq t2 '(tepid hot))]
[(hot) (memq t2 '(warm))]))