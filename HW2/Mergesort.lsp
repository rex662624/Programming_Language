;defun Mymerge (Myarray front mid end)

(defun Mymerge (a low high)
(let 
	;defined variable
	(	
		(b (make-array (+ 1 (- high low))))
		(k -1)
		(mid (floor (+ low high) 2) )
		(tmp)
	)
	;expression
	(setf tmp (+ mid 1))

	; Merging.
	(do ((i low (+ i 0)) (j tmp (+ j 0)))
		((and (> i mid) (> j high)))
		(cond
			((> i mid)
				(progn
					(setf k (+ k 1))
					(setf (aref b k) (aref a j))
					(setf j (+ j 1))
				)
			)
			((> j high)
				(progn
					(setf k (+ k 1))
					(setf (aref b k) (aref a i))
					(setf i (+ i 1))
				)
			)
			((>= (aref a i) (aref a j))
				(progn
					(setf k (+ k 1))
					(setf (aref b k) (aref a j))
					(setf j (+ j 1))
				)
			)
			(t
				(progn
					(setf k (+ k 1))
					(setf (aref b k) (aref a i))
					(setf i (+ i 1))
				)
			)
		)
	)

	(setf k 0)
	(do ((i low (+ i 1))) 
		((> i high))
		(setf (aref a i) (aref b k))
		(setf k (+ k 1))
	)
)

)

(defun mergesort (Myarray front end)
(let	
	(mid)
	(setf mid (floor (+ front end) 2))
	
	(if(< front end)
		(progn
			(mergesort Myarray front mid)
			(mergesort Myarray (+ mid 1) end)
			(MyMerge Myarray front end)	
		)	
	)

)
)

;main function
(let
;--------part1:declare variable-----
	 ((n (read))
	(numbers))
	
;--------part2:expression-----------
    ;set n element array
	(setf numbers(make-array n))

	;store number in array
	(do ((i 0 (+ i 1)))
		((>= i n))                   
		(setf (aref numbers i) (read)) 
    )

	;call merge sort
	(mergesort numbers 0 (- n 1 ))	
	;print sorted array
	(do ((i 0 (+ i 1)))
		((>= i n))
 		(format t "~A " (aref numbers i)) 
	)
	
	(format t "~%")
	
)
