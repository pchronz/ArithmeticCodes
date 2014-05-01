# XXX Reimplement scaleboundaries, avoiding rounding, conversion and BigFloat.
# XXX Implement decodewindow without recursion.
# TODO Use generic functions to be flexible with the source alphabet.
# TODO Single Gaussian probabilistic model.
# TODO tests?
# TODO Decode based on Vector{Uint8, 1} instead of Chars with binary symbols.
# TODO Avod copying where it is not necessary and use subarrays instead.
# TODO Extend function names that manipulate parameters with an exclamation mark
# TODO Re-use calculations. For example in the probabilistc model don't 
#	recompute the whole model if not necessary, instead update it. Use as co-routine?!

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
end
function TwoDIterator(first=false)
	if first
		TwoDIterator(uint(1), uint(1), false, false)
	else
		TwoDIterator(uint(0), uint(0), false, false)
	end
end
abstract DataSource
type TwoDGreyImage <: DataSource
	data::Array{Uint8, 2}
	overflow::Vector{Uint8}
end
function TwoDGreyImage(data::Array{Uint8, 2})
	TwoDGreyImage(data, Array(Uint8, 0))
end
function hor2diterator(src::TwoDGreyImage, iterator::TwoDIterator)
	if iterator.done
		iterator.overflow = true
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
	if !src_it.overflow
		src.data[src_it.h, src_it.w] = sym
	else
		push!(src.overflow, sym)
	end
end

# encode the image pixel by pixel
function encode(src::DataSource, src_it_state::SourceIterator, src_it::Function,
	 calculateboundaries::Function)
	# Encoded vector
	encoded::Vector{Char} = Array(Char, 0)
	# This vector counts the occurences of symbols.
	lower::Uint32 = uint32(0)
	upper::Uint32 = uint32(0) - uint32(1)
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
	lowerend::Uint32 = lower + (uint32(1) << (32 - zeropos))
	assert(lowerend > lower)
	assert(lowerend + (uint(1) << (32 - zeropos)) - 1 <= upper)
	append!(encoded, collect(bits(lowerend)[1:zeropos]))

	encoded
end

function calculateoverlap(lower::Uint32, upper::Uint32)
	# find out which bits are the same
	samebits = ~(lower $ upper)
	regmatch = match(r"(1*)(.*)$", bits(samebits))
	overlaplength = length(regmatch.captures[1])
	overlaplength
end

function decode!(encoded::Vector{Char}, alphabet::Vector{Uint8}, src::DataSource, 
	src_it_state::SourceIterator, src_it::Function, 
	calculateboundaries::Function)
	# Window that holds the open code symbols, until a new source symbol 
		#is identified.
	window::Uint32 = 0
	# Pointer to the current position in the interval window string
	pos = 0
	# Lower and upper boundaries for the symbols.
	lower::Uint32 = 0
	upper::Uint32 = 0 - 1
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
		lower, upper, window, pos = decodewindow!(lower, upper, window, 
			pos, alphabet, src, src_it_state, src_it, 
			calculateboundaries)
	end
	src
end

function decodewindow!(lower::Uint32, upper::Uint32, window::Uint32, pos::Int,
	alphabet::Vector{Uint8}, src::DataSource, src_it_state::SourceIterator, 
	src_it::Function, calculateboundaries::Function)
	# The upper boundary of the window interval.
	# Since we cannot represent "1", we will compare as "<=" according to 
		# our presicion instead of "<".
	windowtop::Uint32 = window + (uint32(1) << (32 - pos)) - uint32(1)
	# Indicator whether the window fit into one of the symbol intervals
	# in a particular iterator over all symbols.
	found::Bool = true
	while found
		found = false
		for sym in [collect(alphabet), EOF]
			low::Uint32
			up::Uint32
			if sym == EOF
				low, up = calculateboundaries(lower, upper, src)
			else
				pushsymbol!(src, src_it_state, uint8(sym))
				low, up = calculateboundaries(lower, upper, src,
					src_it_state)
			end
			# Check whether the symbol interval encloses the window.
			if low <= window && up >= windowtop
				# If it is the EOF symbol, you're done.
				if sym == EOF
					info("Reached EOF symbol.")
					return lower, upper, window, pos
				end
				# Otherwise the symbol is already in the target,
				# and you need to increment the iterator, update
				# the boundaries, and shift both.

				# Increment the iterator.
				src_it_state = src_it(src, src_it_state)

				# Shift the boundaries.
				# Identify the identical bits and shift them out.
				overlaplength = calculateoverlap(low, up)
				# Shift the boundaries to remove the written bits
				low = low << overlaplength
				up = up << overlaplength
				window = window << overlaplength
				pos -= overlaplength
				windowtop = window + (uint32(1) << (32 - pos)) - uint32(1)

				# DEBUG
				push!(lower_dec, low)
				push!(upper_dec, up)

				assert(low < up)
				assert(window >= low)
				assert(window <= up)
				assert(pos >= 0)

				# Update the boundaries.
				lower = low
				upper = up

				# Signal that we need to go even deeper.
				found = true
				break
			end
		end
	end
	lower, upper, window, pos
end

function scaleboundaries(lower::Uint32, upper::Uint32, low::Float64, 
	up::Float64)
	# scale the boundaries
	p::Uint32 = upper - lower
	upper::Uint32 = uint32(lower) + uint32(round(float64(p)*float64(up))) - uint32(1)
	lower::Uint32 = uint32(lower) + uint32(round(float64(p)*float64(low)))
	lower, upper
end

# Static probability model for grayscale images.
function calculateboundariesgraystatic(lower::Uint32, upper::Uint32, 
	src::DataSource, it_state::SourceIterator)
	# Get the most recent symbol.
	sym::Uint8 = src.data[it_state.h, it_state.w]
	
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	p_delta = 1.0/256.0 * (1.0 - p_EOF)
	assert(p_delta*256 == 1.0 - p_EOF)

	# Calculate the floating point boundaries.
	low = sym*p_delta
	up = (sym + 1)*p_delta
	scaleboundaries(lower, upper, low, up)
end
# EOF version
function calculateboundariesgraystatic(lower::Uint32, upper::Uint32, 
	src::DataSource)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	low = 1.0 - p_EOF
	up = 1.0
	scaleboundaries(lower, upper, low, up)
end

# Static probability model.
function calculateboundariesstatic(lower::Uint32, upper::Uint32, 
	src::DataSource, it_state::SourceIterator)
	# Get the most recent symbol.
	sym::Uint8 = src.data[it_state.h, it_state.w]
	
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	p_lo = float64(size(src.data)[1])/float64(num_pixels + 1)
	p_hi = 1 - p_EOF - p_lo

	# Calculate the floating point boundaries.
	if sym == 0
		low = 0.0
		up = p_lo
	elseif sym == 255
		low = p_lo
		up = p_lo + p_hi
	else
		error("The probablistic model cannot cope with symbol $sym.")
	end
	scaleboundaries(lower, upper, low, up)
end
# EOF version
function calculateboundariesstatic(lower::Uint32, upper::Uint32, 
	src::DataSource)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	low = 1.0 - p_EOF
	up = 1.0
	scaleboundaries(lower, upper, low, up)
end

# Laplace's rule of succession as adaptive probabilit model.
function calculateboundarieslaplace(lower::Uint32, upper::Uint32, 
	src::DataSource, it_state::SourceIterator)
	# Get the most recent symbol.
	sym::Uint8
	if !it_state.overflow
		sym = src.data[it_state.h, it_state.w]
	else
		sym = src.overflow[end]
	end
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
	if sym == 0
		low = 0.0
		up = p_lo
	elseif sym == 255
		low = p_lo
		up = 1.0 - p_EOF
	else
		error("The probablistic model cannot cope with symbol $sym.")
	end
	scaleboundaries(lower, upper, low, up)
end
# EOF version
function calculateboundarieslaplace(lower::Uint32, upper::Uint32, 
	src::DataSource)
	num_pixels = prod(size(src.data))
	p_EOF = 1.0/float64(num_pixels + 1)
	low = 1.0 - p_EOF
	up = 1.0
	scaleboundaries(lower, upper, low, up)
end

# Test it.
img_w = 50
img_h = 50

EOF = "EOF"

function testGrayscaleStatic()
	info("Testing the static model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	#img_64::Vector{Int} = rand(DiscreteUniform(0, 255), img_h*img_w)
	img_64::Vector{Int} = rand(Binomial(256, 0.01), img_h*img_w) - 1
	img_8::Vector{Uint8} = map((x)->uint8(x), img_64)
	img = Image(reshape(img_8, img_h, img_w))
	ImageView.display(img)

	tic()
	# Static probabilistic model.
	encoded = encode(TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesgraystatic)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	decoded = decode!(encoded, Uint8[0:255], 
		TwoDGreyImage(zeros(Uint8, size(img.data))), TwoDIterator(true),
		hor2diterator, calculateboundariesgraystatic)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	lower_comp = lower_enc[1:end - (length(lower_enc) - 
		length(lower_dec))] .- lower_dec
	upper_comp = upper_enc[1:end - (length(upper_enc) - 
		length(upper_dec))] .- upper_dec
	assert(sum(lower_comp) + sum(upper_comp) == 0)

	assert(img.data == img_decoded.data)
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

	tic()
	# Static probabilistic model.
	encoded = encode(TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundariesstatic)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Static probabilistic model.
	decoded = decode!(encoded, Uint8[0, 255], 
		   TwoDGreyImage(zeros(Uint8, size(img.data))), 
		   TwoDIterator(true), hor2diterator, calculateboundariesstatic)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
	#info(reduce((x, y)->"$x $y", lower_enc))
	#info(reduce((x, y)->"$x $y", lower_dec))
	#info("\n" * reduce((x, y)->"$x\n$y", map(x->"$(x[1])\t$(x[3])\n$(x[2])\t$(x[4])\n\n", zip(lower_enc, lower_dec, upper_enc, upper_dec))))
	#lower_comp = lower_enc[1:end - (length(lower_enc) - 
	#	length(lower_dec))] .- lower_dec
	#upper_comp = upper_enc[1:end - (length(upper_enc) - 
	#	length(upper_dec))] .- upper_dec
	#assert(sum(lower_comp) + sum(upper_comp) == 0)

	assert(img.data == img_decoded.data)
	info("Test done")
end

function testBinaryLaplace()
	info("Testing the Laplace model on a binary image.")
	# DEBUG
	global lower_enc = Array(Uint32, 0)
	global upper_enc = Array(Uint32, 0)
	global lower_dec = Array(Uint32, 0)
	global upper_dec = Array(Uint32, 0)

	img = Image(ones(Uint8, img_w, img_h).*0)
	img.data[:, 2] = 255
	ImageView.display(img)

	tic()
	# Laplace's rule of succession as adaptive model.
	encoded = encode(TwoDGreyImage(img.data), TwoDIterator(), hor2diterator,
		calculateboundarieslaplace)
	toc()
	info("Encoded length: $(length(encoded))")

	tic()
	# Laplace's rule of succession as adaptive model.
	decoded = decode!(encoded, Uint8[0, 255], 
		TwoDGreyImage(zeros(Uint8, size(img.data))), TwoDIterator(true), 
		hor2diterator, calculateboundarieslaplace)
	toc()
	info("$(img_w*img_h - length(decoded.data)) pixels were missing")
	img_decoded = Image(decoded.data)
	ImageView.display(img_decoded)

	# DEBUG
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
testGrayscaleStatic()

end

