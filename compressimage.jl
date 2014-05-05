# TODO Use binary search on intervals.
# TODO Only compute the necessary values in the interval. When going linearly,
	# compute one by one. When using binary search, compute half by half.
# TODO Make encoding and decoding independent of the type of symbol again.
# TODO Sub-sample recomputation of intervals.
# TODO There are errors on some images. Randomize the sources, save them, and 
	# identify and solve the problem.
# TODO Use generic functions to be flexible with the source alphabet.
# TODO Tests
# TODO Decode based on Vector{Uint8, 1} instead of Chars with binary symbols.
# TODO Avoid copying where it is not necessary and use subarrays instead.
# TODO Extend function names that manipulate parameters with an exclamation mark

module arithmo4astro
using Images
using ImageView
using Distributions

abstract SourceIterator
type TwoDIterator <: SourceIterator
	w::Uint
	h::Uint
	done::Bool
	overflow::Bool
	copy::TwoDIterator
	TwoDIterator(w::Uint, h::Uint, done::Bool, overflow::Bool) = new(w, h, done, overflow)
	TwoDIterator(w::Uint, h::Uint, done::Bool, overflow::Bool, copy::TwoDIterator) = new(w, h, done, overflow, copy)
end
function TwoDIterator(first=false)
	TwoDIterator(uint(0), uint(0), false, false)
end
abstract DataSource
type TwoDGreyImage <: DataSource
	data::Array{Uint8, 2}
	overflow::Vector{Uint8}
end
function TwoDGreyImage(data::Array{Uint8, 2})
	TwoDGreyImage(data, Array(Uint8, 0))
end
function hor2diterator(src::TwoDGreyImage, iterator::TwoDIterator, holdcopy::Bool = false)
	if iterator.done
		iterator.overflow = true
		if holdcopy
			iterator.copy = TwoDIterator(iterator.w, iterator.h, iterator.done, iterator.overflow)
		end
		return iterator
	else
		if iterator.w == 0 && iterator.h == 0
			iterator.w = 1
			iterator.h = 1
			if holdcopy
				iterator.copy = TwoDIterator(iterator.w, iterator.h, iterator.done, iterator.overflow)
			end
			return iterator
		end
		height, width = size(src.data)
		iterator.w += 1
		if iterator.w > width
			iterator.h += 1
			iterator.w = 1
		elseif iterator.w == width && iterator.h == height
			iterator.done = true
		end
		if holdcopy
			iterator.copy = TwoDIterator(iterator.w, iterator.h, iterator.done, iterator.overflow)
		end
		return iterator
	end
end
function pushsymbol!(src::TwoDGreyImage, src_it::TwoDIterator, sym::Uint8)
	if !src_it.overflow
		src.data[src_it.h, src_it.w] = sym
	else
		push!(src.overflow, sym)
	end
end

# encode the image pixel by pixel
function encode(alphabet::Vector{Uint8}, src::DataSource, src_it_state::SourceIterator, src_it::Function,
	 calculateboundaries::Function, refresh_rate::Int)
	# Encoded vector
	encoded::Vector{Char} = Array(Char, 0)
	# This vector counts the occurences of symbols.
	lower::Uint32 = uint32(0)
	upper::Uint32 = uint32(0) - uint32(2)
	lower_fl::Float64
	upper_fl::Float64
	# XXX Make it generic again.
	intervals::Vector{Float64} = Array(Float64, length(alphabet) + 2)
	# Compute the indices for the symbols.
	# XXX Will only work numberic alphabets, unless we first transform the alphabet to a numeric alphabet.
	# XXX What a waste of memory! Consider using a dictionary.
	symbol_idxs::Dict{Uint8, Int} = Dict{Uint8, Int}()
	for (i,s) in enumerate(alphabet)
		symbol_idxs[s] = i
	end
	# A counter for the refresh rate of the probability model.
	refresh_ctr::Int = 0
	while !src_it_state.done
		if refresh_ctr % refresh_rate == 0
			calculateboundaries(intervals, 0, src, src_it_state)
		end
		src_it_state = src_it(src, src_it_state)
		refresh_ctr += 1
		# XXX Make it generic again.
		symbol_idx::Int = symbol_idxs[src.data[src_it_state.h, src_it_state.w]]
		lower_fl = intervals[symbol_idx]
		upper_fl = intervals[symbol_idx + 1]
		@assert lower_fl < upper_fl
		lower, upper = scaleboundaries(lower, upper, lower_fl, upper_fl)
		@assert lower < upper

		overlaplength = calculateoverlap(lower, upper)
		
		# write out the overlapping bits
		append!(encoded, collect(bits(lower)[1:overlaplength]))
		
		# shift the boundaries to remove the written bits
		lower = lower << overlaplength
		upper = upper << overlaplength
		@assert lower < upper

		# DEBUG
		push!(lower_enc, lower)
		push!(upper_enc, upper)
	end
	# Terminate the message
	# Get the interval for the EOF symbol.
	# XXX Make it generic again.
	if refresh_ctr % refresh_rate == 0
		calculateboundaries(intervals, 0, src, src_it_state)
	end
	lower_fl = intervals[end - 1]
	upper_fl = intervals[end]
	lower, upper = scaleboundaries(lower, upper, lower_fl, upper_fl)

	# DEBUG
	push!(lower_enc, lower)
	push!(upper_enc, upper)

	# The largest interval in [lower, upper) is the one that adds 1 bit to 
		# the *third* 0 after the overlap.
	# On the third zero: the first zero on the lower boundary comes directly
		# after the overlapping part because the lower boundary is less
		# than the upper boundary.
	# Since both differ, the lower boundary has to have a 0 there. Then we 
		# need to add a 1 somewhere to be larger than the lower 
		# boundary. This will specify an interval, whose upper boundary 
		# we get by adding another 1 at the same position. Even then 
		# the first zero will remain untouched and thus the upper 
		# boundary of our new interval will be less than the upper 
		# boundary of the symbol boundary. The third zero!
	# XXX What to do if there are no three zeros in the remaining string?
	overlaplength = calculateoverlap(lower, upper)
	# The second 0 after the overlap.
	zeropos = search(bits(lower), '0', overlaplength + 2)
	zeropos = search(bits(lower), '0', zeropos + 1)
	@assert zeropos > overlaplength
	# write out the overlapping bits
	lowerend::Uint32 = lower + (uint32(1) << (32 - zeropos))
	@assert lowerend > lower
	@assert lowerend + (uint(1) << (32 - zeropos)) - 1 <= upper
	append!(encoded, collect(bits(lowerend)[1:zeropos]))

	encoded
end

function calculateoverlap(lower::Uint32, upper::Uint32)
	# find out which bits are the same
	samebits::Uint32 = ~(lower $ upper)
	# TODO Implement binary search
	done::Bool = false
	pos::Int = 1
	while !done && pos <= 32
		res::Uint32 = samebits | (uint32(1) << (32 - pos))
		if res != samebits
			done = true
		else
			pos += 1
		end
	end
	pos - 1
end

function decode!(encoded::Vector{Char}, alphabet::Vector{Uint8}, src::DataSource, 
	src_it_state::SourceIterator, src_it::Function, 
	calculateboundaries::Function, refresh_rate::Int)
	# Window that holds the open code symbols, until a new source symbol 
		#is identified.
	window::Uint32 = 0
	# Pointer to the current position in the interval window string
	pos::Int = 0
	# Lower and upper boundaries for the symbols.
	lower::Uint32 = uint32(0)
	upper::Uint32 = uint32(0) - uint32(2)
	# The counter for the refresh rate for the probability model.
	refresh_ctr::Int = 0
	# Allocate the intervals vector, to allow for binary search.
	intervals::Vector{Float64} = Array(Float64, length(alphabet) + 2)
	# Go through the target symbols
	for tc in encoded
		pos += 1
		if pos > 32
			# TODO I think it should never happen that we enter this
				# here. Prove it!
			error("The decoding window is too small. The position in
				the window is $pos and the window is of type 
				$(typeof(window))")
		end
		# Convert the target symbol to an integer
		t = uint32(tc) - 48
		# Push the codeword symbol into the window
		window = window | (t << (32 - pos))
		# Decode all source symbols using the codeword interval.
		lower, upper, window, pos, refresh_ctr = decodewindow!(intervals, alphabet, lower, upper, window, 
			pos, src, src_it_state, src_it, calculateboundaries, refresh_ctr, refresh_rate)
	end
	src
end

function decodewindow!(intervals::Vector{Float64}, alphabet::Vector{Uint8}, lower::Uint32, upper::Uint32, window::Uint32, pos::Int,
	src::DataSource, src_it_state::SourceIterator, 
	src_it::Function, calculateboundaries::Function, refresh_ctr::Int, refresh_rate::Int)
	# The upper boundary of the window interval.
	# Since we cannot represent "1", we will compare as "<=" according to 
		# our presicion instead of "<".
	windowtop::Uint32 = window + uint32(1 << (32 - pos)) - uint32(1)
	# Floating point versions.
	window_fl::Float64 = float64(window - lower + 1)/float64(upper - lower + 1)
	wintop_fl::Float64 = float64(windowtop - lower + 1)/float64(upper - lower + 1)
	# Indicator whether the window fit into one of the symbol intervals
		# in a particular iterator over all symbols.
	found::Bool = true
	while found
		found = false
		# Management for binary search.
		binsearch_done::Bool = false
		lower_idx::Int = 1
		upper_idx::Int = length(alphabet) + 1
		# Reset the progress indicator
		symbol_progress::Int = 0
		while !binsearch_done
			# Get the next symbol
			idx::Int = int(floor((upper_idx - lower_idx + 1)/2)) + lower_idx
			if refresh_ctr % refresh_rate == 0
				calculateboundaries(intervals, symbol_progress, src, src_it_state, idx)
				symbol_progress = idx
			else
				calculateboundaries(intervals, symbol_progress, src, src_it_state.copy, idx)
				symbol_progress = maximum(idx, symbol_progress)
			end

			low::Uint32
			up::Uint32
			low_fl::Float64 = intervals[idx]
			up_fl::Float64 = intervals[idx + 1]

			@assert low_fl < up_fl
			@assert window_fl < wintop_fl
			# There are five cases to test for binary search:
			# 1 The symbol interval matches the window. --> we're done here.
			# 3 The symbol interval crosses the upper boundary. --> fail.
			# 4 The symbol interval is below the lower boundary. --> Go below.
			# 5 The symbol interval is above the upper boundary. --> Go above.
			# 6 The window encloses the symbol interval. --> fail.
			# Check case 1:
			# Does the interval match the window?
			if low_fl <= window_fl && up_fl >= wintop_fl
				# If it is the EOF symbol, you're done.
				if idx == length(alphabet) + 1
					info("Reached EOF symbol.")
					return lower, upper, window, pos, refresh_ctr
				end
				# Otherwise push the symbol into the target,
				# and you need to increment the iterator, update
				# the boundaries, and shift both.

				# Increment the iterator.
				if refresh_ctr % refresh_rate == 0
					src_it_state = src_it(src, src_it_state, true)
				else
					src_it_state = src_it(src, src_it_state)
				end
				refresh_ctr += 1

				pushsymbol!(src, src_it_state, alphabet[idx])

				low, up = scaleboundaries(lower, upper, low_fl, up_fl)
				@assert window >= low
				@assert windowtop <= up

				# Shift the boundaries.
				# Identify the identical bits and shift them out.
				overlaplength = calculateoverlap(low, up)
				# Shift the boundaries to remove the written bits
				low = low << overlaplength
				up = up << overlaplength
				window = window << overlaplength
				pos -= overlaplength
				windowtop = window + (uint32(1) << (32 - pos)) - uint32(1)

				window_fl = float64(window - low + 1)/float64(up - low + 1)
				wintop_fl = float64(windowtop - low + 1)/float64(up - low + 1)

				# DEBUG
				push!(lower_dec, low)
				push!(upper_dec, up)

				@assert low < up
				@assert window >= low
				@assert window <= up
				@assert pos >= 0

				# Update the boundaries.
				lower = low
				upper = up

				# Signal that we need to go even deeper.
				found = true
				binsearch_done = true
				# Check case 2:
				# Does the interval cross the lower boundary?
			elseif low_fl < window_fl && up_fl >= window_fl
				binsearch_done = true
				# Check case 3:
				# Does the interval cross the upper boundary?
			elseif up_fl > wintop_fl && low_fl <= wintop_fl
				binsearch_done = true
				# Check case 4:
				# Is the interval below the window?
			elseif low_fl < window_fl && up_fl < window_fl
				lower_idx = idx + 1
				# Check case 5:
				# Is the interval above the window?
			elseif low_fl > wintop_fl && up_fl > wintop_fl
				upper_idx = idx - 1
				# Check case 6:
				# Does window enclose the symbol interval?
			elseif low_fl >= window_fl && up_fl <= wintop_fl
				binsearch_done = true
			else
				info("window_fl -> $(window_fl)")
				info("wintop_fl -> $(wintop_fl)")
				info("low_fl -> $(low_fl)")
				info("up_fl -> $(up_fl)")
				error("Reached an impossible state during binary search.")
			end
		end
	end
	lower, upper, window, pos, refresh_ctr
end

function scaleboundaries(lower::Uint32, upper::Uint32, low::Float64, 
	up::Float64)
	# scale the boundaries
	p::Uint32 = upper - lower + 1
	uppr::Uint32 = uint32(lower) + uint32(round((float64(p) + 1.0)*float64(up))) - uint32(1)
	lowr::Uint32 = uint32(lower) + uint32(round((float64(p) + 1.0)*float64(low)))
	lowr, uppr
end

# Static probability model for grayscale images.
# Given the observed values in src up to the iterator compute the probabilites on the alphabet up to the upper idx symbol.
function calculateboundariesgraystatic(intervals::Vector{Float64}, progress::Int, src::DataSource, it_state::SourceIterator, upper_idx::Int=0)
	if upper_idx == 0
		upper_idx = 257
	end
	if progress == upper_idx
		return 
	end

	num_pixels = prod(size(src.data))
	p_EOF = 1.0/(num_pixels + 1)
	p_delta = 1.0/256.0 * (1.0 - p_EOF)
	@assert p_delta*256 == 1.0 - p_EOF

	intervals[1] = 0.0
	for i = (progress + 1):upper_idx
		intervals[i + 1] = intervals[i] + p_delta
	end
	intervals[end] = 1.0
end

# Generalised rule of succession probability model for grayscale images.
# Given the observed values in src up to the iterator compute the probabilites on the alphabet up to the upper idx symbol.
function calculateboundariesgraysuccession(intervals::Vector{Float64}, progress::Int, src::DataSource, it_state::SourceIterator, upper_idx::Int=0)
	if upper_idx == 0 || upper_idx == 257
		upper_idx = 256
	end
	if progress == upper_idx
		return 
	end

	num_pixels = prod(size(src.data))
	p_EOF = 1.0/(num_pixels + 1)

	# Count the number of occurences
	# XXX Performance: Don't re-allocate and recompute all of this each time.
	occurence_ctr::Vector{Uint} = zeros(Uint, 256)
	for h = 1:it_state.h, w = 1:size(src.data, 2)
		occurence_ctr[src.data[h, w] + 1] += 1
		if h == it_state.h && w == it_state.w
			break
		end
	end

	intervals[1] = 0.0
	for i = (progress + 1):upper_idx
		intervals[i + 1] = intervals[i] + (1 - p_EOF)*(occurence_ctr[i] + 1)/(sum(occurence_ctr) + 256)
	end
	intervals[end] = 1.0
end

# Static probability model.
function calculateboundariesstatic(intervals::Vector{Float64}, progress::Int, src::DataSource, it_state::SourceIterator, upper_idx::Int=0)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	p_lo = float64(size(src.data, 1))/float64(num_pixels + 1)
	p_hi = 1 - p_EOF

	intervals[1] = 0.0
	intervals[2] = p_lo
	intervals[3] = p_hi
	intervals[4] = 1.0
end

# Laplace's rule of succession as adaptive probability model.
function calculateboundarieslaplace(intervals::Vector{Float64}, progress::Int, src::DataSource, it_state::SourceIterator, upper_idx::Int=0)
	img_w = size(src.data)[2]

	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)

	# XXX Want this function to run faster? Re-implement using a for-loop.
	# Linearize all elements up to the iterator.
	src_lin::Vector{Uint8} = reshape(src.data[1:it_state.h, :]', 
		int64(it_state.h * img_w), 1)[:, 1]
	# Sum up the occurences up to the iterator.
	num_covered::Int
	if it_state.h == 0
		num_covered = 0
	else
		num_covered = (it_state.h - 1)*img_w + it_state.w
	end
	# Sum and normalize the values.
	# XXX This summation dominates the execution time. Save the state 
		# somewhere and compute the values progressively.
	p_hi = sum(src_lin[1:num_covered]/255)
	p_hi = (p_hi + 1)/(num_covered + 2)*(1.0 - p_EOF)
	p_lo = (1.0 - p_EOF) - p_hi

	intervals[1] = 0.0
	intervals[2] = p_lo
	intervals[3] = 1.0 - p_EOF
	intervals[4] = 1.0
end

# Test it.
img_w = 1000
img_h = 100

EOF = "EOF"

function testgrayscalesuccession()
	info("Testing the succession model on a grayscale image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	#img_64::Vector{Int} = rand(DiscreteUniform(0, 255), img_h*img_w)
	#img_64::Vector{Int} = rand(Binomial(255, 0.5), img_h*img_w)
	img_64::Vector{Int} = 2ones(Int, img_h*img_w)
	img_8::Vector{Uint8} = map((x)->uint8(x), img_64)
	img = Image(reshape(img_8, img_h, img_w))
	#f = open("/Users/pchronz/Desktop/114.txt")
	#img_data = deserialize(f)
	#close(f)
	#img = Image(img)
	ImageView.display(img)

	refresh_rate::Int = 1

	tic()
	# Static probabilistic model.
	encoded = encode(Uint8[0:255], TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesgraysuccession, refresh_rate)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	#Profile.clear()
	decoded = TwoDGreyImage(zeros(Uint8, size(img.data)))
	decode!(encoded, Uint8[0:255], decoded, TwoDIterator(true),
		hor2diterator, calculateboundariesgraysuccession, refresh_rate)
	#Profile.print()
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	@assert sum(lower_comp) + sum(upper_comp) == 0

	@assert img.data == img_decoded.data
	info("Test done")
end

function testGrayscaleStatic()
	info("Testing the static model on a grayscale image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	#img_64::Vector{Int} = rand(DiscreteUniform(0, 255), img_h*img_w)
	#img_64::Vector{Int} = rand(Binomial(255, 0.5), img_h*img_w)
	img_64::Vector{Int} = 2ones(Int, img_h*img_w)
	img_8::Vector{Uint8} = map((x)->uint8(x), img_64)
	img = Image(reshape(img_8, img_h, img_w))
	#f = open("/Users/pchronz/Desktop/114.txt")
	#img_data = deserialize(f)
	#close(f)
	#img = Image(img)
	ImageView.display(img)

	refresh_rate::Int = 1

	tic()
	# Static probabilistic model.
	encoded = encode(Uint8[0:255], TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesgraystatic, refresh_rate)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	#Profile.clear()
	decoded = TwoDGreyImage(zeros(Uint8, size(img.data)))
	decode!(encoded, Uint8[0:255], decoded, TwoDIterator(true),
		hor2diterator, calculateboundariesgraystatic, refresh_rate)
	#Profile.print()
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	@assert sum(lower_comp) + sum(upper_comp) == 0

	@assert img.data == img_decoded.data
	info("Test done")
end

function testBinaryStatic()
	info("Testing the static model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	img = Image(ones(Uint8, img_w, img_h).*255)
	img.data[:, 2] = 0
	ImageView.display(img)

	refresh_rate::Int = 1

	tic()
	# Static probabilistic model.
	encoded = encode(Uint8[0, 255], TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesstatic, refresh_rate)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	decoded = TwoDGreyImage(zeros(Uint8, size(img.data)))
	decode!(encoded, Uint8[0, 255], decoded, TwoDIterator(true),
		hor2diterator, calculateboundariesstatic, refresh_rate)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	#info(reduce((x, y)->"$x $y", lower_enc))
	#info(reduce((x, y)->"$x $y", lower_dec))
	#info("\n" * reduce((x, y)->"$x\n$y", map(x->"$(x[1])\t$(x[3])\n$(x[2])\t$(x[4])\n\n", zip(lower_enc, lower_dec, upper_enc, upper_dec))))
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	@assert sum(lower_comp) + sum(upper_comp) == 0

	@assert img.data == img_decoded.data
	info("Test done")
end

function testBinaryLaplace()
	info("Testing the Laplace model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	img = Image(ones(Uint8, img_w, img_h).*255)
	img.data[:, 2] = 0
	ImageView.display(img)

	refresh_rate::Int = 1

	tic()
	# Laplace's rule of succession as adaptive model.
	encoded = encode(Uint8[0, 255], TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundarieslaplace, refresh_rate)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Laplace's rule of succession as adaptive model.
	decoded = TwoDGreyImage(zeros(Uint8, size(img.data)))
	decode!(encoded, Uint8[0, 255], decoded, TwoDIterator(true),
		hor2diterator, calculateboundarieslaplace, refresh_rate)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	@assert sum(lower_comp) + sum(upper_comp) == 0

	@assert img.data == img_decoded.data
	info("Test done")
end

testBinaryStatic()
testBinaryLaplace()
testGrayscaleStatic()
testgrayscalesuccession()

end

