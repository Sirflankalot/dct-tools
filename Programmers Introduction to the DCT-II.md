# A Programmers Introduction to the DCT-II

###What is a DCT-II?

The discrete cosine transform is one of the most important algorithms in modern encoding.
Used in formats such as .JPG, .MP4, and even the upcoming Daala codec, it's prevalence makes it a useful transform to know.
What is does is convert the actual values within a NxN square of the image into a bunch of cosine waves that when recombined together will reform the original image.
The advantage of this is that you can then get rid of the extremely high frequency sign waves easily and without changing too much of the image.

###Rationale

When I was doing research to teach myself how to implement a DCT I found a couple of different articles and implementations.
The problem was that no matter which article I followed my code could never worked and produced the results they provided.
This was most likely because I didn't fully understand what they were asking me to do and they weren't very programmer friendly articles.
Therefore, I figured it would be a good idea to take all I have learned trying to implement the DCT-II and share it in a more step by step, programmer oriented article.
In this article I'm going to provide example code in both C++ and Python to increase the amount of people that will be able to use this, as well as some tools to make developing your own implementation of DCT-II easier. 
I will also include visualizations as the amount of nesting that is needed can be confusing without being able to see things. 

##The Math

While the DCT-II has a fairly daunting formula behind it, once you break it down into it's component parts, it is actually fairly straight forward and simple to implement.
The actual formula for a DCT-II is as follows:

<p align="center">![DCT-II Formula](https://wikimedia.org/api/rest_v1/media/math/render/svg/dce6d60796ea026a5a7564418d130effde90d9cf)

Even for someone who likes math, this formula is not very easy to parse.

This isn't even the 2D implementation of DCT-II which is:

<p align="center">![2D DCT-II Formula](https://wikimedia.org/api/rest_v1/media/math/render/svg/4a3639270e74e1b69ec593b56e752474e0bf365f)

This is even worse.
Despite how daunting this looks, you can actually cheat by using two sets of one dimensional DCT-IIs.