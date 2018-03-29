(defun mergesort (numbers)
 (return-from mergesort numbers))




;main function
(let
	;part1:declare variable
	 ((n (read))
      (numbers))
	;part2:expression
 (setf numbers
  (do ( (i 0 (+ i 1)) (tmp nil))
	((>= i n)	(reverse tmp))
       (setf tmp (cons (read) tmp))
  )
 )

    (format t "~{~A ~}~%" (mergesort numbers)))
