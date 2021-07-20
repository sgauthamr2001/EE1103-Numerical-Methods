# Hamming Distance 

The [Hamming distance](https://www.geeksforgeeks.org/hamming-distance-two-strings/) is one easy way of matching two bit streams. We use this assignment to understand the complexity of a quantum erasure channel, and with a QBER as high as 25%. Start by : 

Using arguments (inputs) to main()  
Generating random numbers, estimate the Hamming distance between arrays of random bits, without errors  
For your submission, your program will have 3 inputs  
- N: the total number of bits ~ 1 million (more the better, but start small)
- P: the fraction of bits that reach the receiver ~ 0.1
- Q: the fraction of bits that have an error ~ 0.25

The assignment has four parts  
- Setup A and B arrays, create C as a subset of B.
- A will use C (with unknown offset and some % errors), plotting Hamming distance versus offset. Whenever the Hamming distance is a minimum, you have the true offset.
- Show that, in the presence of errors, your ability to determine the location of min(Hdist) improves as length(B') increases.
- Check the randomness of the secret key that you have simulated.
