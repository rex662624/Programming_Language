(defun Mymerge (Myarray front mid end)
(let 

	;defined variable
	(	
		(left (make-array (+ 2 (- mid front))));front~mid
		(right (make-array(+ 1 (- end mid))));mid+1~end
		(leftmax (+ 1 (- mid front)))
		(rightmax (+ 1 (- end (+ 1 mid))))
		(rightidx 0)
		(leftidx 0)
		(i 0)
		(j 0)
	)

	;expression
	
	;把array[front]~array[mid]放進 LeftSub[]
	(
	do((i front (+ i 1)) (j 0 (+ 1 j)))
	((> i mid))
	 (setf (aref  left j) (aref Myarray i)) 
;	(format t "~A " (aref left j)) 
	)
	(setf (aref  left leftmax) -1)	

	;把array[mid+1]~array[end]放進 RightSub[]
	(
 	do((i  (+ 1 mid) (+ i 1)) (j 0 (+ 1 j)))
 	((> i end))
	(setf (aref  right j) (aref Myarray i))
;	(format t "~A " (aref right j)) 
	)
	(setf (aref  right rightmax) -1)

	(do ((i front (+ i 1)))
	((> i end))
		(cond
		(
	(and (not (= (aref left leftidx) -1)) (or (= (aref right rightidx) -1) (<= (aref left leftidx) (aref right rightidx))))
			  	(progn
				(setf (aref Myarray i) (aref left leftidx))
				(setf leftidx (+ leftidx 1))	 
				)	
		)			
			(t ;else
				(setf (aref Myarray i) (aref right rightidx))
				(setf rightidx (+ rightidx 1))
			)
		) 
	)
)
)

(defun mergesort (Myarray front end)
(let	
	(
		(mid (floor (+ front end) 2))
	)
	
	(if(< front end)
		(progn
			(mergesort Myarray front mid)
			(mergesort Myarray (+ mid 1) end)
			(MyMerge Myarray front mid end)	
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
