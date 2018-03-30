(defun fib1(n)
  (cond
    ((eq n 1) 0)
    ((eq n 2) 1)
    ((+ (fib1(- n 1)) (fib1(- n 2))))
	)
)

(defun fib2(n a result)
    (cond
    ((eq n 1) a)
    ((eq n 2) result)
    ((fib2 (- n 1) result (+ a result)))
	)

)

; MAIN FUNCTION
(let 
	(
	(n (read))
	)
	
	(trace fib1)
	(trace fib2)
	(fib1 n)	
	;(format t "orig: ~A ~%" (fib1 n))
	(format t "-------------------------~%" )
	(fib2 n 0 1)
	;(format t "tail: ~A ~%" (fib2 n 0 1))
	
)

