# TODO Grayscale images.
# TODO Use generic functions to be flexible with the source alphabet.
# TODO Single Gaussian probabilistic model.
# TODO tests?
# TODO Decode based on Vector{Uint8, 1} instead of Chars with binary symbols.
# TODO Avod copying where it is not necessary and use subarrays instead.
# TODO Extend function names that manipulate parameters with an exclamation mark
# TODO Re-use calculations. For example in the probabilistc model don't 
#	recompute the whole model if not necessary, instead update it. Use as co-routine?!
# XXX Why doesn't the algorithm work properly for 2x2?

module arithmo4astro
using Images
using ImageView

abstract SourceIterator
type TwoDIterator <: SourceIterator
	w::Uint
	h::Uint
	done::Bool
end
function TwoDIterator(first=false)
	if first
		TwoDIterator(uint(1), uint(1), false)
	else
		TwoDIterator(uint(0), uint(0), false)
	end
end
abstract DataSource
type TwoDGreyImage <: DataSource
	data::Array{Uint8, 2}
end
function hor2diterator(src::TwoDGreyImage, iterator::TwoDIterator)
	if iterator.done
		return iterator
	else
		if iterator.w == 0 && iterator.h == 0
			iterator.w = 1
			iterator.h = 1
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
		return iterator
	end
end
function pushsymbol!(src::TwoDGreyImage, src_it::TwoDIterator, sym::Uint8)
	src.data[src_it.h, src_it.w] = sym
end

# encode the image pixel by pixel
function encode(src::DataSource, src_it_state::SourceIterator, src_it::Function,
	 calculateboundaries::Function)
	# Encoded vector
	encoded::Vector{Char} = Array(Char, 0)
	# This vector counts the occurences of symbols.
	lower::Uint64 = 0
	upper::Uint64 = 0 - 1
	while !src_it_state.done
		src_it_state = src_it(src, src_it_state)
		lower, upper = calculateboundaries(lower, upper, src, 
			src_it_state)
		assert(lower < upper)

		overlaplength = calculateoverlap(lower, upper)
		
		# write out the overlapping bits
		append!(encoded, collect(bits(lower)[1:overlaplength]))
		
		# shift the boundaries to remove the written bits
		lower = lower << overlaplength
		upper = upper << overlaplength
		assert(lower < upper)

		# DEBUG
		push!(lower_enc, lower)
		push!(upper_enc, upper)
	end
	# Terminate the message
	# Get the interval for the EOF symbol.
	lower, upper = calculateboundaries(lower, upper, src)

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
	assert(zeropos > overlaplength)
	# write out the overlapping bits
	lowerend = lower + (uint64(1) << (64 - zeropos))
	assert(lowerend > lower)
	assert(lowerend + (uint(1) << (64 - zeropos)) - 1 <= upper)
	append!(encoded, collect(bits(lowerend)[1:zeropos]))
	
	encoded
end

function calculateoverlap(lower::Uint64, upper::Uint64)
	# find out which bits are the same
	samebits = ~(lower $ upper)
	regmatch = match(r"(1*)(.*)$", bits(samebits))
	overlaplength = length(regmatch.captures[1])
	overlaplength
end

function decode!(encoded::Vector{Char}, src::DataSource, 
	src_it_state::SourceIterator, src_it::Function, 
	calculateboundaries::Function)
	# Window that holds the open code symbols, until a new source symbol 
		#is identified.
	window::Uint64 = 0
	# Pointer to the current position in the interval window string
	pos = 0
	# Lower and upper boundaries for the symbols.
	lower::Uint64 = 0
	upper::Uint64 = 0 - 1
	# Go through the target symbols
	for tc in encoded
		pos += 1
		if pos > 64
			# TODO I think it should never happen that we enter this
				# here. Prove it!
			error("The decoding window is too small. The position in
				the window is $pos and the window is of type 
				$(typeof(window))")
		end
		# Convert the target symbol to an integer
		t = uint64(tc) - 48
		# Push the codeword symbol into the window
		window = window | (t << (64 - pos))
		# Decode all source symbols using the codeword interval.
		lower, upper, window, pos = decodewindow!(lower, upper, window, 
			pos, src, src_it_state, src_it, calculateboundaries)
	end
	src
end

function decodewindow!(lower::Uint64, upper::Uint64, window::Uint64, pos::Int,
	src::DataSource, src_it_state::SourceIterator, src_it::Function,
	calculateboundaries::Function)
	# The upper boundary of the window interval.
	# Since we cannot represent "1", we will compare as "<=" according to 
		# our presicion instead of "<".
	windowtop = window + (uint64(1) << (64 - pos)) - 1
	for sym in [locol, hicol, EOF]
		if sym == EOF
			low, up = calculateboundaries(lower, upper, src)
		else
			pushsymbol!(src, src_it_state, uint8(sym))
			low, up = calculateboundaries(lower, upper, src, 
				src_it_state)
		end
		# Test whether our window is within the symbol's interval
		if low <= window && up >= windowtop
			if sym == EOF
				info("Reached EOF symbol.")
				return lower, upper, window, pos
			end
			src_it_state = src_it(src, src_it_state)
			# Identify the identical bits and shift them out.
			overlaplength = calculateoverlap(low, up)
			
			# shift the boundaries to remove the written bits
			low = low << overlaplength
			up = up << overlaplength
			window = window << overlaplength
			pos -= overlaplength
			assert(low < up)
			assert(window >= low)
			assert(window <= up)
			assert(pos >= 0)

			# DEBUG
			push!(lower_dec, low)
			push!(upper_dec, up)

			return decodewindow!(low, up, window, pos, src, 
				src_it_state, src_it, calculateboundaries)
		end
	end
	lower, upper, window, pos
end

function scaleboundaries(lower::Uint64, upper::Uint64, low::Float64, 
	up::Float64)
	# scale the boundaries
	p = upper - lower
	upper = lower + uint64(round(BigFloat(p)*BigFloat(up))) - 1
	lower = lower + uint64(round(BigFloat(p)*BigFloat(low)))
	lower, upper
end

# Static probability model.
function calculateboundariesstatic(lower::Uint64, upper::Uint64, 
	src::DataSource, it_state::SourceIterator)
	# Get the most recent symbol.
	sym::Uint8 = src.data[it_state.h, it_state.w]
	
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	p_lo = float64(size(src.data)[1])/float64(num_pixels + 1)
	p_hi = 1 - p_EOF - p_lo

	# Calculate the floating point boundaries.
	if sym == locol
		low = 0.0
		up = p_lo
	elseif sym == hicol
		low = p_lo
		up = p_lo + p_hi
	else
		error("The probablistic model cannot cope with symbol $sym.")
	end
	scaleboundaries(lower, upper, low, up)
end
# EOF version
function calculateboundariesstatic(lower::Uint64, upper::Uint64, 
	src::DataSource)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	low = 1.0 - p_EOF
	up = 1.0
	scaleboundaries(lower, upper, low, up)
end

# Laplace's rule of succession as adaptive probabilit model.
function calculateboundarieslaplace(lower::Uint64, upper::Uint64, 
	src::DataSource, it_state::SourceIterator)
	# Get the most recent symbol.
	sym::Uint8 = src.data[it_state.h, it_state.w]
	img_w = size(src.data)[2]

	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)

	# Linearize all elements up to the iterator.
	src_lin::Vector{Uint8} = reshape(src.data[1:it_state.h, :]', 
		int64(it_state.h * img_w), 1)[:, 1]
	# Sum up the occurences up to the iterator.
	num_covered = (it_state.h - 1)*img_w + it_state.w - 1
	# Sum and normalize the values.
	# XXX This summation dominates the execution time. Save the state 
		# somewhere and compute the values progressively.
	p_hi = sum(src_lin[1:num_covered]/255)
	p_hi = (p_hi + 1)/(num_covered + 2)*(1.0 - p_EOF)
	p_lo = (1.0 - p_EOF) - p_hi

	# Calculate the floating point boundaries.
	if sym == locol
		low = 0.0
		up = p_lo
	elseif sym == hicol
		low = p_lo
		up = 1.0 - p_EOF
	else
		error("The probablistic model cannot cope with symbol $sym.")
	end
	scaleboundaries(lower, upper, low, up)
end
# EOF version
function calculateboundarieslaplace(lower::Uint64, upper::Uint64, 
	src::DataSource)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	low = 1.0 - p_EOF
	up = 1.0
	scaleboundaries(lower, upper, low, up)
end

# Test it.
img_w = 100
img_h = 100

locol = 0
hicol = 255
EOF = "EOF"

function testBinaryStatic()
	info("Testing the static model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint64, 0)
	global upper_enc = Array(Uint64, 0)
	global lower_dec = Array(Uint64, 0)
	global upper_dec = Array(Uint64, 0)

	img = Image(ones(Uint8, img_w, img_h).*hicol)
	img.data[:, 2] = locol
	ImageView.display(img)

	tic()
	# Static probabilistic model.
	encoded = encode(TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesstatic)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	decoded = decode!(encoded, TwoDGreyImage(zeros(Uint8, size(img.data))), 
		TwoDIterator(true), hor2diterator, calculateboundariesstatic)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	#info("length(lower_enc) -> $(length(lower_enc))")
	#info("length(upper_enc) -> $(length(upper_enc))")
	#info("length(lower_dec) -> $(length(lower_dec))")
	#info("length(upper_dec) -> $(length(upper_dec))")
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	assert(sum(lower_comp) + sum(upper_comp) == 0)

	assert(img.data == img_decoded.data)
	info("Test done")
end

function testBinaryLaplace()
	info("Testing the Laplace model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint64, 0)
	global upper_enc = Array(Uint64, 0)
	global lower_dec = Array(Uint64, 0)
	global upper_dec = Array(Uint64, 0)

	img = Image(ones(Uint8, img_w, img_h).*hicol)
	img.data[:, 2] = locol
	ImageView.display(img)

	tic()
	# Laplace's rule of succession as adaptive model.
	encoded = encode(TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundarieslaplace)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Laplace's rule of succession as adaptive model.
	decoded = decode!(encoded, TwoDGreyImage(zeros(Uint8, size(img.data))), 
		TwoDIterator(true), hor2diterator, calculateboundarieslaplace)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	#info("length(lower_enc) -> $(length(lower_enc))")
	#info("length(upper_enc) -> $(length(upper_enc))")
	#info("length(lower_dec) -> $(length(lower_dec))")
	#info("length(upper_dec) -> $(length(upper_dec))")
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	assert(sum(lower_comp) + sum(upper_comp) == 0)

	assert(img.data == img_decoded.data)
	info("Test done")
end

testBinaryStatic()
testBinaryLaplace()

end

