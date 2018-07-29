R Inference and Modeling Course
================

Parameters and Estimates
------------------------

### Sampling Models

Let's create model for the following scenario:

There is an urn filled with beads, some red and some blue. The challenge is to guess the spread between the two colors. To mimic the expense of running polls, each draw costs 10 cents; and to mimic the competition between media companies vying for your attention, the reward for an accurate poll is 25 dollars.

Here are the rules:

1.  If you provide an interval that contains the true proportion, you get half of what you paid to sample and move on to the second part of the competition.

2.  The entry with the smallest interval is named the winner.

The dslabs package has a function for random draws from an urn.

``` r
take_poll(25)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-2-1.png) Note that we have described a simple model for opinion polls (blue beads represents Democrats and red beads are Republicans). We want to predict the proportion of blue beads in the urn.

Some definitions:

> Population: all of the beads in the urn (from which we will sample)

> Parameter: the proportion of beads in the urn of a given color

> Sample: the beads we checked

The task of statistical inference is predict the parameter, *p*, using the observed data in the sample. Note that the proportion of blue beads is *p*, while for red beads it is 1 − *p*. The spread is |*p* − (1 − *p*)| or 2*p* − 1.

> Estimate: summary of observed data that is informative about the parameter of interest

### The Sample Average

Each draw is a random variable *X*. For *n* draws, each draw is represented by *X*<sub>1</sub>, ..., *X*<sub>*n*</sub>.

The average of these random variables is $\\bar{X} = \\frac{X\_1,...,X\_n}{n}$.

We know that the expected value of the sum of *n* draws is *n* times the average of the values in the urn, which must be the proportion *p*. Since we can't find out *p* directly, we will try to estimate it.

### Polling Versus Forecasting

If a poll is conducted four months before an election, it is estimating *p* for that moment, not for election day. Forecasters try to build tools that model how opinions vary accross time, and try to predict the election day result from early polls.

### Properties of Our Estimate

When we multiply $\\bar{X}$ by *n*, we get a sum of independent draws, so we know that our probability rules apply. We can see that the expected value is the following:

$$E(n\\bar{X}) = n \\times p$$
 Dividing by *n*, we can see that the expected value is *p*.

We can compute the standard error with the same method we used in the probability course:

$$SE\[X\] = |a-b|\\sqrt{p(1-p)}$$
 Because we are dividing by the sum *n* for our average of the independent draws, we add *n* to our equation. *a* and *b* are 1 and zero, respectively, se we can simplify.

$$SE\[X\] = |a-b|\\sqrt{p(1-p)/n} = \\sqrt{p(1-p)/n}$$
 Thus we know that the expected value of our poll is *E*\[*X*\]=*p*, and that the standard error will decrease proportional to the square root of *n*. This tells us that with a large enough poll, our estimate converges to *p*.

If we take a large enough poll for our standard error to be, say, 0.01, we can be quite sure of the winner. But how large a poll do you need for that narrow of a range?

We don't actually know *p*, so we can't compute the standard error. For illustrative purposes, let's compute it assuming *p* = 0.51.

``` r
p = 0.51
Ns <- seq(1,20000)
SE <- sqrt(p*(1-p)/Ns)
data.frame(sample_size = Ns, standard_error = SE) %>% ggplot(aes(sample_size, standard_error)) + geom_line() + scale_x_continuous(trans = "log10") + geom_vline(aes(xintercept = 10000))
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-3-1.png)

We would need about 10,000 people to get the standard error down to 0.01. Usually, polls range in sample size from 500 to 3500 people. A poll with a 1000 people where the true proportion is 0.51 has a standard error of 1.5%

Interestingly, the standard error is largest around *p* = 0.5 for a given *N*.

``` r
N <- 25
p <- seq(0,1,length.out = 100)
se <- sqrt(p*(1-p)/N)
plot(p,se)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-4-1.png) With the following graphs, you can see how sample size reduces the SE while the pattern along proportions stays the same.

``` r
for(sample_sizes in c(25,100,1000)){
  se <- sqrt(p*(1-p)/sample_sizes)
  plot(p, se, ylim = c(0,max(se)))
}
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-5-1.png)![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-5-2.png)![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-5-3.png)

It's important to note that the standard error is largest for *p* = 0.5.

The CLT in Practice
-------------------

### Estimating the Standard Error

The Central Limit Theorem tells us that the distribution of a sum of draws is approximately normal.

*X*<sub>1</sub> + *X*<sub>2</sub> + ... + *X*<sub>*N*</sub>
 We also know that dividing by a random variable by a non-random constant will still produce a normal distribution.

$$\\frac{X}{a} \\approx N(\\frac{\\mu}{a},\\frac{\\sigma}{a})$$

This implies that the average of the draws, $\\bar{X}$, is approximately normal.

> Importantly, we know that the expected value of $\\bar{X}$ is *p*, and that the standard error of $\\bar{X}$ is $\\sqrt{p(1-p)/N}$.

Suppose we want to know the probability that we are 1% point from point (i.e. that we made a very good estimate).

In mathematical termrs, our question amounts to the following: $Pr(|\\bar{X}-p) \\le 0.01$.

Recalling what we know about intervals on a normal curve, we can see that the probability above can be expressed as:

$$Pr(\\bar{X}) \\le p + 0.01) - Pr(\\bar{X}) \\le p - 0.01) $$

Let's use the same trick we learned in the previous module: subtracting the expected value and dividing by the stndard error for every term:

$$Pr(\\bar{X}) \\le p + 0.01) - Pr(\\bar{X}) \\le p - 0.01)$$
 becomes:

$$Pr \\bigg(\\frac{\\bar{X} - E\[\\bar{X}\]}{SE\[\\bar{X}\]} \\le \\frac{(p + 0.01) - E\[\\bar{X}\]}{SE\[\\bar{X}\]}\\bigg) - Pr\\bigg(\\frac{\\bar{X} - E\[\\bar{X}\]}{SE\[\\bar{X}\]} \\le \\frac{(p - 0.01) - E\[\\bar{X}\]}{SE\[\\bar{X}\]}\\bigg)$$
 This gives us the standard normal variable, which we call *Z*, on the left side of each propability.

$$Pr \\bigg(Z \\le \\frac{(p + 0.01) - E\[\\bar{X}\]}{SE\[\\bar{X}\]}\\bigg) - Pr\\bigg(Z \\le \\frac{(p - 0.01) - E\[\\bar{X}\]}{SE\[\\bar{X}\]}\\bigg)$$
 From earlier, we know that $E\[\\bar{X}\] = p$, and that $SE\[\\bar{X}\] = \\sqrt{p(1-p)/N}$. Let's substitute in those terms above. Note that *p* becomes $E\[\\bar{X}\]$ which then cancels.

$$Pr \\bigg(Z \\le \\frac{0.01}{\\sqrt{p(1-p)/N}}\\bigg) - Pr\\bigg(Z \\le \\frac{-0.01}{\\sqrt{p(1-p)/N}}\\bigg)$$

Now, can we compute the probablity? Not yet. In the equation above, we don't know *p*. The CLT tells us that we can use an estimate of the standard error that replaces the true proportion *p* in the above equation with $\\bar{X}$. We denote that this is an estimate of the standard with a hat over SE:

$$\\hat{SE} = \\sqrt{\\bar{X}(1-\\bar{X})/N}$$
 Plugging the above into our previous equation, we get:

$$Pr \\bigg(Z \\le \\frac{0.01}{\\sqrt{\\bar{X}(1-\\bar{X})/N}}\\bigg) - Pr\\bigg(Z \\le \\frac{-0.01}{\\sqrt{\\bar{X}(1-\\bar{X})/N}}\\bigg)$$
 We call the above a **plug-in estimate**.

Now, we can compute the probabilities using the available data. Let's say a random draw of beads (our earlier example) had 12 blue beads and 13 red beads. That's a proportion of red beads of 0.48.

``` r
x_hat <- 12/25
se <- sqrt(x_hat*(1-x_hat)/25)
se
```

    ## [1] 0.09991997

Now we can compute the probability being within 1% of the true proportion, *p*, using `pnorm()`.

``` r
pnorm(0.01/se) - pnorm(-0.01/se)
```

    ## [1] 0.07971926

We can see that the probability of our estimate being within 1% of the true proportion is about 8%. Thus, it is unlikely that we will be within our desired range of the proportion. Later, we'll learn how to determine a sample size large enough to make informative estimates and probabilities.

### Margin of Error

The **margin of error** is two times the standard error, which we learned to estimate above.

$$MoE = 2\\hat{SE}$$

The margin of error gives the probability that $\\bar{X}$ is within two standard errors from *p*, which can be represented as follows:

*P**r*(*Z* ≤ 2)−*P**r*(*Z* ≤ −2)

We know from our normal curve that the probability that our estimate is within two standard errors of *p* is about 95%.

``` r
pnorm(2) - pnorm(-2)
```

    ## [1] 0.9544997

It's worth noting that the choice of two standard errors (and thus a 95% probability) is somewhat arbitrary. In our example case, the margin of error is $MoE = 2 \\times \\hat{SE} \\approx 0.2$. Clearly, a margin of error of 20% from our sample of 25 isn't particularly useful.

To provide potentially useful estimates, political polls tend to have sample sizes in the between 700 and 3,500. We can see how this gives a more powerful estimate with a hypothetical scenario where 2,000 beads are drawn and the $\\bar{X}$ is still 0.48.

``` r
x_hat <- 0.48
se <- sqrt(x_hat*(1-x_hat)/2000)
2*se
```

    ## [1] 0.02234278

We can see that the margin of error for our hypothetical is 2%, so we could feel more confident that there are fewer blue beads than red ones. Of course, we might not get the same {X} if we actually ran the simulation.

### Simulating our CLT

Running a Monte Carlo simulation as we have before is complicated by the fact that we don't know *p*. However, we can pick values of *p* and use them to show that the theory used above works well in practice.

Let's pick *p* = 0.45. Let's simulate polling 1000 beads (or people) 10,000 times.

``` r
p <- 0.45
N <- 1000
B <- 10000
X_hat <- replicate(B, {
  X <- sample(c(0,1), size=N, replace = TRUE, prob = c(1-p, p))
  mean(X)
})
mean(X_hat)
```

    ## [1] 0.4498756

``` r
sd(X_hat)
```

    ## [1] 0.0156117

The theory tells us that $\\bar{X}$ has an approximately normal distribution with expected value 0.45 and standard error ~1.5%.

Our simulation confirms this. If we take the means of our $\\hat{X}$s above, we get 0.45. The standard deviation of our simulated $\\hat{X}$ gives the same standard error as well.

A histogram and a Q-Q plot of the simulation show that the assumption of normality is fair as well.

``` r
library(gridExtra)
p1 <- data.frame(X_hat = X_hat) %>% ggplot(aes(X_hat)) +
  geom_histogram(binwidth = 0.005, color = "black")
p2 <- data.frame(X_hat = X_hat) %>% ggplot(aes(sample = X_hat)) +
  stat_qq(dparams = list(mean = mean(X_hat), sd = sd(X_hat))) +
  geom_abline() +
  ylab("X_hat") +
  xlab("Theoretical Normal")
grid.arrange(p1,p2, nrow = 1)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-11-1.png) You would see similar results for different values of *N* and *p*!

### The Spread

Our competition above is to predict the spread, not the proportion. Since there are only two possibilities, the spread is simply *p* − (1 − *p*) or 2*p* − 1. Substituting $\\bar{X}$ for *p* we can see that our estimate of the spread will be $2\\bar{X}-1$, and that our $SE\[2\\bar{X}\]$ will be $2\\hat{SE}\[\\bar{X}\]$, since subtracting 1 has no effect on the standard error. In our first example with *N* = 25, $\\bar{X} = 0.48$ and $\\hat{SE}\[\\bar{X}\] \\approx 0.2$. For the spread, $|2\\bar{X}-1| = 0.04$, and $2\\hat{SE} \\approx 2\*0.2 \\approx 0.4$. Clearly, a 40% interval for the spread is not useful for practical purposes, thus our use of much larger sample sizes for most cases (e.g. political polling).

The standard error of the spread can be expressed as the square root of the average squared distance $(\\bar{X} - p)^2$. In practice, that would be the `sqrt(mean(errors)^2)`, where `errors` represents $\\bar{X}-p$.

### Bias in Polling

If we were to poll 100,000 people in order to find a binary proportion (as we have been), our theoretical standard error would max out at about 0.3% -- a very accurate poll! So why don't we see more polls with this many people?

First off, a poll this large would be very expensive. But, perhaps more importantly, polling is much more complicated than the theory suggests. Polling is usually done by phone, which can leave out people who do not have phones. We also don't actually know who is going to vote, so even if our margin of error is very small; it may not be exactly right that our expected value is *p*. We call this phenomenon **bias**.

CIs and p-Values
----------------

### Confidence Intervals

We can use confidence intervals to find the interval that has a given likelihood of containing *p* -- often, we use an inteval witha 95% likelihood of containing *p*.

We can see CIs in action with `geom_smooth()`.

``` r
data("nhtemp")
data.frame(year = as.numeric(time(nhtemp)), temperature=as.numeric(nhtemp)) %>%
  ggplot(aes(year, temperature)) +
  geom_point() +
  geom_smooth() +
  ggtitle("Average Yearly Temperatures in New Haven")
```

    ## `geom_smooth()` using method = 'loess'

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-12-1.png)

A confidence interval is usually defined as 2 standand errors around the expected value, $\\bar{X}$. We can represent that as follows:

$$CI = \[\\bar{X} - 2\\hat{SE}(\\bar{X}),\\space \\bar{X} +2\\hat{SE}(\\bar{X})\]$$

We want to know the probability that the above interval contains the true proportion, *p*. It's important to note that both values in the interval are random variables, and as such vary randomly each time they are simulated. We can calculate the probability that *p* is within the interval above as follows:

$$ Pr(\\bar{X} - 2\\hat{SE}(\\bar{X}) \\le p \\le  \\bar{X} +2\\hat{SE}(\\bar{X}))$$
 Subtracting $E\[\\bar{X}\]$ (recall that $E\[\\bar{X}\] = p$) and dividing by $\\hat{SE}\[\\bar{X}\]$, this becomes:

$$ Pr\\bigg(-2 \\le \\frac{\\bar{X}-p}{\\hat{SE}(\\bar{X})} \\le  2\\bigg) $$
 The middle term is our standard normal variable, *Z*. We can calculate probability that *Z* falls into the interval above using `pnorm()`.

``` r
pnorm(2) - pnorm(-2)
```

    ## [1] 0.9544997

As you can see, there is an approximately 95% chance that the true proportion, *p*, falls into the CI. If we wanted to compute a 99% CI instead, we can determine with `qnorm()` and `pnorm()` which *z* satisfies the following equation:
*P**r*(−*z* ≤ *Z* ≤ *z*)=0.99

``` r
z <- qnorm(0.995)
pnorm(qnorm(0.995))
```

    ## [1] 0.995

``` r
pnorm(1-qnorm(0.995))
```

    ## [1] 0.05753257

``` r
pnorm(z) - pnorm(-z)
```

    ## [1] 0.99

We can use this approach for any percentile, *q*, using the following equation:

1 − (1 − *q*)/2
 Note that above we used `qnorm(0.995)` to find the appropriate *z*. This is because `qnorm()` calculates the z value over the interval −inf to *z*. To find an interval covering 99% of values, we want to exclude the 1% of values farthest from $E\[\\bar{X}\]$ using the equation above. Entering the desired percentile above will return the value to be entered into `qnorm()` to find the appropriate *z*.

### Simulating CIs

Let's run a Monte Carlo simulation to confirm that *p* is actually contained within our theoretical 95% confidence interval 95% of the time.

``` r
p <- 0.45
N <- 1000
B <- 10000
inside <- replicate(B, {
  X <- sample(c(0,1), N, replace = TRUE, prob = c(1-p, p))
  X_hat <- mean(X)
  SE_hat <- sqrt(X_hat*(1-X_hat)/N)
  between(p, X_hat - 2*SE_hat, X_hat + 2*SE_hat)
})
mean(inside)
```

    ## [1] 0.9548

Here, you can see that in 0.95% of cases, *p* is inside the interval that is generated using 1000 random samples.

``` r
sample_CI <-replicate(100, {
  X <- sample(c(0,1), 1000, replace = TRUE, prob = c(1-p, p))
  X_hat <- mean(X)
  SE_hat <- sqrt(X_hat*(1-X_hat)/N)
  c(X_hat, X_hat - 2*SE_hat, X_hat + 2*SE_hat)
})
sample_CI <- t(sample_CI)
CI_data <- data.frame(mean = sample_CI[,1], low = sample_CI[,2], high = sample_CI[,3])
CI_data <- CI_data %>% mutate(inside = p >= low & p <= high)
```

    ## Warning: package 'bindrcpp' was built under R version 3.2.5

``` r
CI_data %>% ggplot() + geom_segment(aes(x = low, xend=high,y=1:100, yend = 1:100, color = inside), size = 1, show.legend = FALSE) + geom_vline(aes(xintercept = 0.45), size = 1.1, color = "darkgrey") + geom_point(aes(x = mean, y = 1:100)) + xlab("Simulated Interval") + ylab("Polls")+ ggtitle("100 Simulated CIs") + theme_minimal()
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-16-1.png)

The graph above shows a sampling of a hundred simulated intervals, where most contain *p*, but a few do not (and marked in red).

### The Correct Language

It is important to remember that it is the intervals that are random. not *p*. Thus, the 95% refers to the probability that the confidence interval falls on top of *p*. It is technically incorrect to say that there is a 95% chance that *p* falls within a certain interval.

### Power

In the earlier case where we sampled 25 beads to estimate a proportion, we saw that our CI of the spread included zero. We can think of power as the probability of detecting a spread different from 0, which increases as we increase our sample size and shrink the corresponding standard error.

### p-Values

Suppose we take a random sample of 100 beads and draw 52 blue ones. The spread, 2*p* − 1 is 4%. However, we know that it's possible that the true spread is zero, and this result is the result of chance.

> The *null hypothesis* is the skeptic's hypothesis that there is in fact no difference.

In this case, that would mean that the spread is zero.

The p-value tells us likely it is to see a 4% spread given that there true spread is zero.

> The p-value asks: How likely is it to see a value this large when the null hypothesis is true?

In this case, we can ask that question by asking: what's the probability that the absolute value of $\\bar{X} - 0.5$ is greater than 0.02.

$$Pr(|\\bar{X} - 0.5|) &gt; 0.02$$
 Under the null hypothesis, we know that the following expression ($\\bar{X}$ minus p and divided by the standard error) is a standard normal variable.

$$\\sqrt{N}\\frac{\\bar{X}- 0.5}{0.5(1-0.5)}$$
 With this, we can derive an equation for the p-value.

$$Pr\\bigg(\\sqrt{N}\\frac{\\bar{X}- 0.5}{0.5(1-0.5)} &gt; \\sqrt{N}\\frac{0.02}{0.5(1-0.5)}\\bigg)$$

$$Pr\\bigg(Z &gt; \\sqrt{N}\\frac{0.02}{0.5(1-0.5)}\\bigg)$$

``` r
N = 100
z = sqrt(N)*0.02/0.5
1 - (pnorm(z)-pnorm(-z))
```

    ## [1] 0.6891565

In this case, the p-value is 0.689, meaning that there is a 68.9% chance that you would see 52 beads on 100 draws if the true proportion was 0.5. This is very likely! So, we are not convinced that the spread is nonzero.

If a 95% confidence interval does not include 0, we know that the p-value must be smaller than 0.05. It's important to note that the p-value reports a probability, but says nothing about the significance of the result in the context of the problem.

Statistical Models
------------------

### Poll Aggregators

Poll aggregators use the results of multiple polls to increase the precision with which they predict outcomes. Let's take an illustrative example wherein we plot the confidence intervals of 12 simulated polls (say, of Obama or Romney):

``` r
d <- 0.039
Ns <- c(1298, 533, 1342, 897, 774, 254, 812, 324, 1291, 1056, 2172, 516)
p <- (d + 1)/2
confidence_intervals <- sapply(Ns, function(N) {
  X <- sample(c(0,1), size = N, replace = TRUE, prob = c(1-p, p))
  X_hat <- mean(X)
  SE_hat <- sqrt(X_hat*(1-X_hat)/N)
  2*c(X_hat, X_hat - 2*SE_hat, X_hat + 2*SE_hat)-1
})
polls <- data.frame(poll=1:ncol(confidence_intervals), 
                    t(confidence_intervals),
                      sample_size = Ns)
names(polls) <- c("poll", "estimate", "low", "high", "sample_size")
polls %>% ggplot() + geom_segment(aes(x = low, xend=high,y=poll, yend = poll), size = 1, show.legend = FALSE) + geom_vline(aes(xintercept = 0.039), size = 1.1, color = "darkgrey") +xlab("Spread (true proportion marked with grey line)")
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-18-1.png) Nearly all of the polls above include the actual spread, 3.9%, but nearly all also include 0. Thus, most pollsters would have to say that the election is a tossup. However, we can aggregate these results to get a more precise sense of the spread.

We don't have access to the raw poll data, but we can use our knowledge of distributions to combine them into a theoretical poll a sample size equal to the sum of that of the individual polls.

We can construct an estimate of the spread, or *d*, with an average weighted for the sample sizes.

``` r
d_hat <- polls %>% 
  summarize(avg = sum(estimate*sample_size)/sum(sample_size)) %>% .$avg #recall that estimate refers to the estimated spread
d_hat
```

    ## [1] 0.03079244

We can use our estimate of the spread, $\\hat{d}$, to estimate the proportion voting for Obama, *p*, and the margin of error.

``` r
p_hat <- (1+d_hat)/2
moe <- 2*1.96*sqrt(p_hat*(1-p_hat)/sum(polls$sample_size)) #multiply SE by 1.96 to get the components of our interval, and then by 2 to get the full margin of error

round(d_hat*100,1)
```

    ## [1] 3.1

``` r
round(moe*100,1)
```

    ## [1] 1.8

``` r
polls
```

    ##    poll    estimate          low       high sample_size
    ## 1     1 0.003081664 -0.052430810 0.05859414        1298
    ## 2     2 0.039399625 -0.047162727 0.12596198         533
    ## 3     3 0.050670641 -0.003854336 0.10519562        1342
    ## 4     4 0.059085842 -0.007575547 0.12574723         897
    ## 5     5 0.025839793 -0.046024718 0.09770430         774
    ## 6     6 0.015748031 -0.109727568 0.14122363         254
    ## 7     7 0.066502463 -0.003528404 0.13653333         812
    ## 8     8 0.006172840 -0.104936155 0.11728183         324
    ## 9     9 0.010069713 -0.045590498 0.06572993        1291
    ## 10   10 0.058712121 -0.002727455 0.12015170        1056
    ## 11   11 0.018416206 -0.024490623 0.06132304        2172
    ## 12   12 0.011627907 -0.076411231 0.09966705         516

``` r
polls %>% rbind(c(13, d_hat, d_hat - moe/2, d_hat + moe/2, sum(polls$sample_size))) %>% mutate(agg = sample_size < 10000) %>% ggplot() + geom_segment(aes(x = low, xend=high,y=poll, yend = poll, color = agg), size = 1, show.legend = FALSE) + geom_vline(aes(xintercept = 0.039), size = 1.1, color = "darkgrey") +xlab("Spread (true proportion marked with grey line)")
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-20-1.png) , we can then predict that the spread will be 3.1% with a margin of error of 1.8%. Along with being more precise, this makes us much more confident that Obama is ahead. The aggregate result above is marked in red.

Note that this is just a simulation, and that the actual forecasting is more complicated, involving complex statistical modeling.

### Poll Data and Pollster Bias

Here's a simplified version of a model for the popular vote that is reminiscent of the methods used at FiveThirtyEight for the 2016 presidential election.

``` r
polls_us <- polls_us_election_2016 %>% filter(state == "U.S." & enddate >= "2016-10-31" & (grade %in% c("A+", "A", "A-","B") | is.na(grade))) %>% 
  mutate(spread = rawpoll_clinton/100 - rawpoll_trump/100) 
```

Let's assume that there are only two partes, so that the proportions for Clinton and Trump are *p* and 1 − *p*, respectively. We'll call the spread between the two proportions *d*.

The theory we learned tells us that the *d*<sub>*p**o**l**l*</sub> of each poll is a random variable with a probability distribution that is approximately normal. *E*\[*d*<sub>*p**o**l**l*</sub>\] is the election night spread, *d*, and *S**E*\[*d*<sub>*p**o**l**l*</sub>\] is $2\\sqrt{p(1-p)/N}$.

``` r
d_hat_2 <- polls_us %>% summarize(d_hat_2 = sum(spread * samplesize) / sum(samplesize)) %>% .$d_hat
p_hat_2 <- (d_hat+1)/2
moe_2 <- 2*1.96*sqrt(p_hat_2*(1-p_hat_2)/sum(polls_us$samplesize))
d_hat_2
```

    ## [1] 0.0262393

``` r
moe_2
```

    ## [1] 0.004115113

Using this data, we get a spread of 3.1% with a margin of error of 0.4%. On election night, we find out that the actual percentage was 2.1% -- outside our interval. Was this just bad luck?

A look at the histogram of the spreads from our set of polls reveals that the spreads are not normally distributed and that the standard error is larger than we estimated. The theory is not quite working here.

``` r
polls_us %>%
  ggplot(aes(spread)) +
  geom_histogram(color = "black", binwidth = 0.01)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-23-1.png)

Here's a table showing how many polls each pollster conducted in the last week.

``` r
polls_us %>% group_by(pollster) %>% summarize(n())
```

    ## # A tibble: 16 x 2
    ##    pollster                                                   `n()`
    ##    <fct>                                                      <int>
    ##  1 ABC News/Washington Post                                       7
    ##  2 Angus Reid Global                                              1
    ##  3 CBS News/New York Times                                        2
    ##  4 Fox News/Anderson Robbins Research/Shaw & Company Research     2
    ##  5 Google Consumer Surveys                                        2
    ##  6 IBD/TIPP                                                       8
    ##  7 Insights West                                                  1
    ##  8 Ipsos                                                          6
    ##  9 Marist College                                                 1
    ## 10 Monmouth University                                            1
    ## 11 Morning Consult                                                1
    ## 12 NBC News/Wall Street Journal                                   1
    ## 13 Selzer & Company                                               1
    ## 14 The Times-Picayune/Lucid                                       8
    ## 15 USC Dornsife/LA Times                                          8
    ## 16 YouGov                                                         3

We can see that some are polling much more frequently than others. Let's filter for pollsters that poll at least 6 times and plot their spreads.

``` r
polls_us %>% group_by(pollster) %>%
  filter(n() >= 6) %>%
  ggplot(aes(pollster, spread)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
polls_us %>% group_by(pollster) %>% 
  summarize(se = 2*sqrt(p_hat_2*(1-p_hat_2)/median(samplesize)))
```

    ## # A tibble: 16 x 2
    ##    pollster                                                        se
    ##    <fct>                                                        <dbl>
    ##  1 ABC News/Washington Post                                   0.0265 
    ##  2 Angus Reid Global                                          0.0295 
    ##  3 CBS News/New York Times                                    0.0291 
    ##  4 Fox News/Anderson Robbins Research/Shaw & Company Research 0.0288 
    ##  5 Google Consumer Surveys                                    0.00627
    ##  6 IBD/TIPP                                                   0.0333 
    ##  7 Insights West                                              0.0326 
    ##  8 Ipsos                                                      0.0225 
    ##  9 Marist College                                             0.0326 
    ## 10 Monmouth University                                        0.0365 
    ## 11 Morning Consult                                            0.0260 
    ## 12 NBC News/Wall Street Journal                               0.0279 
    ## 13 Selzer & Company                                           0.0354 
    ## 14 The Times-Picayune/Lucid                                   0.0196 
    ## 15 USC Dornsife/LA Times                                      0.0183 
    ## 16 YouGov                                                     0.0165

We can see that these polls have standard errors that we would expect, but their expected values differ. The theory tells us that all of the polls should have the same expected value, so this tells us that the different pollsters are subject to "house effects", or bias.

We'll use data-driven models rather independent urn-draw model to improve the accuracy of our predictions.

### Data-Driven Models

For each pollster, let's collect their last reported result before the election.

``` r
one_poll_per_pollster <- polls_us %>% group_by(pollster) %>% 
  filter(enddate == max(enddate)) %>% ungroup()
one_poll_per_pollster %>%
  ggplot(aes(spread)) + geom_histogram(binwidth = 0.01)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-27-1.png) We will use a new model to more accurately estimate the spread. We assume that the expected value of our spread is the actual spread, 2*p* − 1. Our new urn contains polls from all possible pollsters rather than zeros and ones, so it now contains continuous numbers from -1 to 1. This mean that the *standard deviation* of the urn is no longer $\\sqrt{(}p(1-p)}$. Our new standard error for our average includes pollster-to-pollster variability as well as the sample variability. Regardless, the standard deviation is now an *unknown parameter* represented by *σ*.

> The sample average is a random variable with expected value *μ* and standard error $\\frac{\\sigma}{\\sqrt{N}}$.

The central limit theorem still applies to the average of the drawn values (which are independent variables). For a large enough N, the probability distribution is approximately normal with expected value *d* and standard deviation *σ*. If we consider N=15 large enough, we can use this to construct a confidence interval.

We can estimate *σ* (the unobserved standard deviation) using the *sample standard deviation*, defined as follows:

$$s = \\frac{1}{N - 1}\\sqrt{\\sum\_{i=1}^N(X\_i-\\bar{X}^2)}$$
 The `sd()` function in R computes the sample standard deviation using the above formula, so we can compute for our data as follows:

``` r
sd(one_poll_per_pollster$spread)
```

    ## [1] 0.02357651

We use the CLT to compute a confidence interval with the following code:

``` r
results <- one_poll_per_pollster %>% 
  summarize(avg = mean(spread), se = sd(spread)/sqrt(length(spread))) %>% 
  mutate(start = avg - 1.96*se, end = avg + 1.96*se)
round(results*100,1)
```

    ##   avg  se start end
    ## 1 2.9 0.6   1.8 4.1

Our new interval still strongly suggests that Clinton is ahead in the popular vote, while consider pollster-to-pollster variability. Can we estimate her probability at this point? Not yet, since *d* is a fixed parameter. To do that, we need Bayesian statistics.

Bayesian Statistics
-------------------

In the previous sections, *p* has been a fixed parameter, so it doesn't make sense to estimate its probability. With Bayesian statistics, we can treat it as a random variable. We'll also discuss multilevel modeling approaches such as *hierarchical modeling* wherein different types of variability (e.g. sampling, pollster-to-pollster, week-to-week, election-to-election) are considered separately.

### Bayes' Theorem

Let's use a cystic fibrosis test to demonstrate the theorem. If a test is 99% accurate, we can say that the chance the test is positive given that one has the disease is 0.99:

*P*(+|*D* = 1)=0.99

The converse is also true:

*P*(−|*D* = 0)=0.99
 If we select a random person and they test positive, what are the chances they have the disease given that the cystic fibrosis rate is 1:3900, or:

*P**r*(*D* = 1)≈0.000256
 In effect, we are looking for answer to the following expression, which asks what the probability is that D = 1 given that the test is positive:

$$P(D = 1 | +) = \\space ?$$
 To answer this, we use *Bayes' Theorem*:

$$Pr(A|B) = \\frac{Pr(A \\space and \\space B)}{Pr(B)}$$
 Using the multiplication rule, the numerator is split:

$$Pr(A|B) = \\frac{Pr(B|A)Pr(A)}{Pr(B)}$$
 Let's apply this to our cystic fibrosis example:

$$PP(D = 1 | +) = \\frac{Pr(+ | D = 1) \\times Pr(D = 1)}{Pr(+)}$$
Expanding the denominator to include only terms we know:
$$PP(D = 1 | +) = \\frac{Pr(+ | D = 1) \\times Pr(D = 1)}{Pr(+|D = 1) \\times Pr(D = 1) + Pr(+|D=0)Pr(D=0)}$$

Now, that we have all values we know, we plug in our values:

$$PP(D = 1 | +) = \\frac{0.99 \\times 0.000256}{0.99 \\times 0.000256 + 0.01 \\times 0.000975} \\approx 0.02$$
 Despite having a 99% accuracy, the probability that the patient has the disease given that the test is positive is only 2%!

This is because we have to account for the rarity of the disease in the population.

Let's do a Monte Carlo simulation to demonstrate this:

First, we sample from 100,000 "people" whether they have the disease using the prevalence of cystic fibrosis.

``` r
prev <- 0.00025 #the prevalence of the disease in the population
N <- 100000 #our sample population (getting the test)
outcome <- sample(c("Disease","Healthy"), N, replace = TRUE, 
                  prob = c(prev,1-prev))
N_D <-sum(outcome == "Disease")
N_D #only a few "people" have the disease since it is so rare
```

    ## [1] 19

``` r
N_H <-sum(outcome == "Healthy")
N_H #the vast majority of the sampled "patients" are healthy
```

    ## [1] 99981

Let's create a vector with hypothetical test results between the healthy and CF groups given 99% accuracy of the test. Here, we sample for a correct or incorrect test with the appropriate probability. We can then generate a table of outcomes.

``` r
accuracy <- 0.99 #the accuracy of the test
test <- vector("character",N) #create a vector of length N with mode = "character"
#sample to generate test outcomes
test[outcome=="Disease"]  <- sample(c("+","-"), N_D, replace=TRUE, 
                                    prob = c(accuracy, 1 - accuracy)) 
test[outcome=="Healthy"]  <- sample(c("-","+"), N_H, replace=TRUE, 
                                    prob = c(accuracy, 1 - accuracy))
cf <- data.frame(outcome = outcome, test = test) #not necessary for table (for #s in description)
table(outcome, test)
```

    ##          test
    ## outcome       -     +
    ##   Disease     0    19
    ##   Healthy 99027   954

We can see that all 19 people with the disease tested positive, but so did a small proportion but large number of healthy people. The number of people with CF testing positive (19) divided by the total testing positive (973) provides our expected proportion of about 2%.

### Bayes' Theorem in Practice

Let's use a baseball example to demonstrate the usefulness of Bayes' theorem in practice. Jose Iglesias is a professional baseball player who had a 0.450 batting average after his first 20 at-bats in 2013.

With the techniques we have learned up to know (*frequentist statistics*), the best we can do is provide a confidence interval.

We treat the average as the true proportion, *p*, and compute the standard error:

$$SE = \\sqrt{\\frac{0.450(1-0.450)}{20}} = 0.111$$
 We calculate our confidence interval as 0.228 to 0.678 with a mean of 0.450. Immediately, we can see that this isn't particularly useful for two reasons:

1.  the interval is very large, ranging from far better than the record to below average.

2.  a 0.450 batting average at the end of the season would be world record!

It turns out that the average batting avg. over the 2010, 2011, and 2012 seasons among all MLB players was 0.275 with a standard deviation of 0.027 -- much less then Iglesias' average after 20 at-bats. So is Jose lucky or the best batter in 50 years? Bayesian statistics will help us decipher that, and better predict how he will do down the line.

### The Hierarchical Model

We use a model with two levels to predict how Iglesias will perform. First, each player is assigned a natural ability to hit, *p*. We assume that *p* is normal and (based on the data above) has an expected value of 0.275 and a standard error 0.027. The second level of variability accounts for luck when batting.

At each at-bat, a player a chance *p* of getting a hit. The CLT tells us that the observed average, *Y*, has a normal distribution where *E*\[*Y*\]=*p* and $SE\[Y\] = \\sqrt(p(1-p)/N)$ where N is the number of at-bats.

We write this as follows:

*p* *N*(*μ*, *τ*<sup>2</sup>) describes randomness in picking a player *Y*|*p* *N*(*p*, *σ*<sup>2</sup>) describes randomness in performance of a particular player

with *μ* = 0.275, *τ* = 0.027, and *σ*<sup>2</sup> = *p*(1 − *p*)/*N*.

> Note that there are two levels: 1) player to player variability, and 2) variability of batting luck by the individual. The first level is the *prior distribution* and the second is the *sampling distribution*

Let's use this model to predict *p*, Jose's innate ability, or "true" batting average.

*p* *N*(0.275, 0.27<sup>2</sup>) (note that *μ* is the average of all players) *Y*|*p* *N*(*p*, 0.111<sup>2</sup>) (note that sigma is the standard error of the observed average for Jose)

The continuous version of Bayes' rule can be used to derive the *posterior probability function*, which is the probability *p* assuming we observe *Y* = *y*. The expected value in this case is:

*E*(*p*|*Y* = *y*)=*B**μ* + (1 − *B*)*y* = *μ* + (1 − *B*)(*y* − *μ*)
 where $B = \\frac{\\sigma^2}{\\sigma^2 +\\tau^2}$

Applying this in Jose's case,

*E*(*p*|*Y* = *y*)=0.275 + (1 − *B*)(0.450 − 0.275)
$$B = \\frac{0.111^2}{0.111^2 +0.027^2}$$
*E*(*p*|*Y* = 0.450)≈0.285

The standard error can be shown to be

$$SE(p|y)^2 = \\frac{1}{1/\\sigma^2 + 1/\\tau^2} = \\frac{1}{1/0.111^2 + 1/0.027^2} \\approx 0.00069 $$
 The standard deviation is the square root of the above, or 0.026.

This is an empirical Bayesion approach, since we used data to find the prior probability. We can now report a *95% credible interval* by reporting an interval centered at the mean with a 95% chance of occuring. In our case, this turns out to be 0.285 ± 0.052 (2 times the standard deviation).

Note that way B is calculated makes intuitive sense. If the standard deviation of the observed average is much larger than the standard deviation of the prior probability (*σ* &gt; &gt;*τ*) and thus very imprecise, then B = 1. When B = 1, the expected value of *p* given *Y* ignores *Y* entirely and has expected value *μ* (the prior expected value). If we don't trust our observed data, we go with the historical result. In the opposite case when *τ* &gt; &gt;*σ* and B = 0, we ignore the past and rely on our observed result.

Election Forecasting
--------------------

### Bayesian Approach

Pollsters attempt to forecast who will an election by providing a probability, e.g. Obama has a 91% chance of winning. Let's use a Bayesian model to predict how likely Clinton is to win the popular vote, where Clinton winning is parameter *d*.

We write it like this:

*d* *N*(*μ*, *τ*) describes our best guess had we not seen any polling data

$\\bar{X} | d ~ N(d,\\sigma)$ describes randomness due to sampling and the pollster effect

A popular approach uses so-called fundamentals such as the state of the economy to compute *μ*, but we will just call it 0 (which just says we don't know who will win). We'll set *τ* to 0.035 since the average spread has been about 3.5%.

``` r
mu <- 0
tau <- 0.035
sigma <- results$se
Y <- results$avg
B <- sigma^2/(sigma^2 + tau^2)
posterior_mean <- B*mu + (1-B)*Y
posterior_se <- sqrt(1/(1/sigma^2 + 1/tau^2))
posterior_mean
```

    ## [1] 0.02832788

``` r
posterior_se #standard error
```

    ## [1] 0.005812285

``` r
posterior_mean + c(-1.96, 1.96)*posterior_se #the 95% credible interval
```

    ## [1] 0.01693580 0.03971996

Above, we calculate a 95% credible interval using the standard error to get a range where the spread is likely to fall. It's more useful, though, to calculate how likely the spread is greater than zero (i.e. CLinton's chances). We do that with `pnorm()`.

``` r
1 - pnorm(0, posterior_mean, posterior_se)
```

    ## [1] 0.9999995

This tells us that Clinton's chances are greater than 99.9%! What accounts for the difference between this FiveThirtyEight's estimate of 81.4%? Our model does not account for the *general bias* of pollsters which can swing from a couple percentage points in one direction to the other between elections. We can't know which way it swings until the election has happened, but we can include a term to account for this additional variability.

### Mathematical Representations of Models

Let's model a pollster without general bias. He collects several polls of sample size *N*, and observes measurements of spread *X*<sub>1</sub>, ..., *X*<sub>*j*</sub> These random variables have expected value d and standard error $2\\sqrt(p(1-p)/N)$.

*X*<sub>*j*</sub> = *d* + *ϵ*<sub>*j*</sub>
 The index *J* represents the different polls, while *ϵ*<sub>*j*</sub> explains poll-to-poll variability due to sampling error. We assume its average is zero, and its standard error is $2\\sqrt(p(1-p)/N)$.

If *d* is 2.1 and N = 2000, we can simulat *J* = 6 as follows:

``` r
set.seed(3)
J<- 6
N <- 2000
d <- .021
p <- (d + 1)/2
X <- d + rnorm(J,0,2*sqrt(p*(1-p)/N))
X
```

    ## [1] -0.0005047417  0.0144603685  0.0267854043 -0.0047567709  0.0253768717
    ## [6]  0.0216734433

Now let's expand to consider multiple pollsters. We add another index *I* for each pollster (5 in our example).

We use *X*<sub>*i**j*</sub> with *i* representing the pollster and *j* representing the j-th pollster. Here we choose 5 pollsters each with 6 polls.

*X*<sub>*i**j*</sub> = *d* + *ϵ*<sub>*i**j*</sub>
 To simulate this, we loop through the pollsters using the `sapply()` function:

``` r
d <- .021
I <- 5
J <- 6
N <- 2000
X <- sapply(1:I, function(i){ ##give sapply vector, open function
  d + rnorm(J,0,2*sqrt(p*(1-p)/N)) ##define function for sapply to apply to vector
})
df_2 <- data.frame(spreads = c(sapply(1:5, function(i){X[,i]})), pollsters = c(sapply(1:5, function(i){rep(i, 6)})))
ggplot() + geom_point(aes(df_2$pollsters, df_2$spreads)) + xlab("Pollsters") + ylab("Spreads") +ylim(c(-0.06, 0.08))
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-35-1.png) This doesn't capture the actual data, which shows variability in the means of the different pollsters. We have to add an additional random variable to account for that. We will use *h*<sub>*i*</sub> to represent the house effect of the i-th pollster.

*X*<sub>*i*, *j*</sub> = *d* + *h*<sub>*i*</sub> + *ϵ*<sub>*i*, *j*</sub>
 Assuming *σ*<sub>*h*</sub> is 0.025, we simulate the data as follows:

``` r
I <- 5
J <- 6
N <- 2000
d <- .021
p <- (d + 1)/2
h <- rnorm(I, 0, 0.025)
X2 <- sapply(1:I, function(i){
  d + h[i] + rnorm(J,0,2*sqrt(p*(1-p)/N))
})
df_3 <- data.frame(spreads = c(sapply(1:5, function(i){X2[,i]})), pollsters = c(sapply(1:5, function(i){rep(i, 6)})))
ggplot() + geom_point(aes(df_3$pollsters, df_3$spreads)) + xlab("Pollsters") + ylab("Spreads") +ylim(c(-0.06, 0.08))
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-36-1.png)

Note that while this looks more like the actual data, where *h* affects each pollsters' mean, we assume that the mean of the house effects is zero. But this is often not the case in the actual data. Thus, we have to add another random variable *b*, which we'll give *σ*<sub>*b*</sub> = 0.025 from the historical data:

*X*<sub>*i*, *j*</sub> = *d* + *b* + *h*<sub>*i*</sub> + *ϵ*<sub>*i*, *j*</sub>

Further note that this implies the following:

$$\\bar{X} = d + b + \\frac{1}{N}\\sum\_{i = 1}^N X\_i$$
 Which has standard deviation:

$$\\sqrt{\\sigma^2/N + \\sigma\_b^2}$$
 &gt;Imporantly, since *b* is the same in every measurement, increasing N does not reduce its variance. No matter how many polls you take, this bias is not reduced.

If we redo our calculation above considering *b*, we get a result much closer to FiveThirtyEight's 81%:

``` r
mu <- 0
tau <- 0.035
sigma <- sqrt(results$se^2 + .025^2)
Y <- results$avg
B <- sigma^2 / (sigma^2 + tau^2)

posterior_mean <- B*mu + (1-B)*Y
posterior_se <- sqrt( 1/ (1/sigma^2 + 1/tau^2))

1 - pnorm(0, posterior_mean, posterior_se)
```

    ## [1] 0.8197348

### Predicting the Electoral College

We've focused on predicting the popular vote, but as we know, the presidential election is not decided by a majority vote, but by the Electoral College.

Here are the top 5 states ranked by electoral votes:

``` r
results_us_election_2016 %>% top_n(5, electoral_votes)
```

    ##          state electoral_votes clinton trump others
    ## 1   California              55    61.7  31.6    6.7
    ## 2      Florida              29    47.8  49.0    3.2
    ## 3     Illinois              20    55.8  38.8    5.4
    ## 4     New York              29    59.0  36.5    4.5
    ## 5 Pennsylvania              20    47.9  48.6    3.6
    ## 6        Texas              38    43.2  52.2    4.5

Let's use a series of filter steps to select high quality statewide polls from the week before the election.

``` r
results <- polls_us_election_2016 %>%
  filter(state!="U.S." & 
           !grepl("CD", state) & #returns a logical based on the pattern 
           enddate >="2016-10-31" & 
           (grade %in% c("A+","A","A-","B+") | is.na(grade))) %>%
  mutate(spread = rawpoll_clinton/100 - rawpoll_trump/100) %>%
  group_by(state) %>%
  summarize(avg = mean(spread), sd = sd(spread), n = n()) %>%
  mutate(state = as.character(state))
```

We can see from a table of the results that Florida, North Carolina, Ohio, etc. have average spreads very close to zero. These are the battleground states.

``` r
results %>% arrange(abs(avg))
```

    ## # A tibble: 47 x 4
    ##    state               avg     sd     n
    ##    <chr>             <dbl>  <dbl> <int>
    ##  1 Florida         0.00356 0.0163     7
    ##  2 North Carolina -0.00730 0.0306     9
    ##  3 Ohio           -0.0104  0.0252     6
    ##  4 Nevada          0.0169  0.0441     7
    ##  5 Iowa           -0.0197  0.0437     3
    ##  6 Michigan        0.0209  0.0203     6
    ##  7 Arizona        -0.0326  0.0270     9
    ##  8 Pennsylvania    0.0353  0.0116     9
    ##  9 New Mexico      0.0389  0.0226     6
    ## 10 Georgia        -0.0448  0.0238     4
    ## # ... with 37 more rows

We can add the number of electoral votes for each state with `left_join()`.

``` r
results <- left_join(results, results_us_election_2016, by = "state")
results_us_election_2016 %>% filter(!state %in% results$state)
```

    ##                  state electoral_votes clinton trump others
    ## 1               Alaska               3    36.6  51.3   12.2
    ## 2         Rhode Island               4    54.4  38.9    6.7
    ## 3              Wyoming               3    21.9  68.2   10.0
    ## 4 District of Columbia               3    90.9   4.1    5.0

``` r
results <- results %>%
  mutate(sd = ifelse(is.na(sd), median(results$sd, na.rm=TRUE), sd)) #replaces the sd with the     median sd where there is only one poll
```

Note that the above 3 states and DC do not appear in `results$state` because no polls were taken there, since the results are pretty much known (Rhode Island and DC to the Democrats and the other two to the Republicans).

Let's create a Monte Carlo simulation to assess the probability that Clinton wins the electoral college.

First, we assume for simplicity's sake that we don't have any previous information about the mean (prior spread of zero) and just take the average variance of historical results (2% or sd of 0.02). Using the equations we learned, we calculate the posterior mean and standard error.

``` r
mu <- 0
tau <- 0.02
results %>% mutate(sigma = sd/sqrt(n), 
                   B = sigma^2 / (sd^2 + tau^2),
                   posterior_mean = B*mu + (1-B)*avg,
                   posterior_se = sqrt( 1/ (1/sigma^2 + 1/tau^2))) %>%
  arrange(abs(posterior_mean)) %>%
  ggplot(aes(avg, posterior_mean, size = n)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0)
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-42-1.png) As you can see from the plot above that the means produced from a larger number of polls are less influenced by the push towards zero that results from using an estimate based on a posterior. A small poll showing a difference of 0.2 in the "avg", for example, is shifted toward 0 to ~0.1 in the posterior mean. The larger polls (which tend to have smaller spreads) don't shift as much (i.e. stay on the identity line).

Now we repeat this 10,000 times using `rnorm()` for our Monte Carlo simulation to estimate how many electoral votes Clinton is likely to win. We add 7 at the end to account for Rhode Island and DC, for which we have no polls, but are sure the Democrats will win.

``` r
mu <- 0
tau <- 0.02
clinton_EV <- replicate(1000, { 
  results %>% mutate(sigma = sd/sqrt(n), 
                   B = sigma^2 / (sigma^2 + tau^2),
                   posterior_mean = B*mu + (1-B)*avg,
                   posterior_se = sqrt( 1/ (1/sigma^2 + 1/tau^2)),
                   simulated_result = rnorm(length(posterior_mean), posterior_mean, posterior_se), #creates a vector of distributed means for each state excluding RI, AK, WI, and DC
                   clinton = ifelse(simulated_result>0, electoral_votes, 0)) %>% #if 
    summarize(clinton = sum(clinton)) %>% #if the spread > 0, then Clinton gets all the electoral votes for that state, otherwise she gets zero; her total is summed
    .$clinton + 7## 7 for Rhode Island and D.C. (Alaska and WI are assumed to be losses, i.e. 0)
})
mean(clinton_EV>269) ##in what proportion of the simulated results does Clinton have enough EC votes to win?
```

    ## [1] 0.998

This model gives Clinton an over 99% chance of winning. The Princeton Election Commission made a similar estimate, which was clearly quite off. This model ignores the general bias of polls during a given election, which was between 1-2% in 2016. Pollsters became over-confident because of the large numbers of polls taken in battleground states which made standard errors appear small.

``` r
data_frame(c_EV = clinton_EV) %>% ggplot(aes(c_EV)) + geom_histogram(binwidth = 1) + geom_vline(aes(xintercept = 269))
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-44-1.png) Let's recompute after taking into account the general bias.

``` r
tau <- 0.02
bias_sd <- 0.03 #we assume that the general bias contributes an additional 3% variance (1.5% up or down)
clinton_EV_2 <- replicate(1000, {
  results %>% mutate(sigma = sqrt(sd^2/n  + bias_sd^2),  ##note the bias_sd term which is not diminished by poll number!
                   B = sigma^2 / (sigma^2 + tau^2),
                   posterior_mean = B*mu + (1-B)*avg,
                   posterior_se = sqrt( 1/ (1/sigma^2 + 1/tau^2)),
                   simulated_result = rnorm(length(posterior_mean), posterior_mean, posterior_se),
                   clinton = ifelse(simulated_result>0, electoral_votes, 0)) %>% 
    summarize(clinton = sum(clinton) + 7) %>% .$clinton ## 7 for Rhode Island and D.C.
})
mean(clinton_EV_2>269)
```

    ## [1] 0.837

With this additional adjustment, we can see that the likelihood has declined substantially, and that the variability of simulated outcomes has increased.

``` r
data_frame(c_EV = clinton_EV_2) %>% ggplot(aes(c_EV)) + geom_histogram(binwidth = 1) + geom_vline(aes(xintercept = 269))
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-46-1.png)

### Forecasting

How important are polls taken several weeks before the election?

To make sure the variability we observe is not due to pollster effects, let’s study data from one pollster:

``` r
one_pollster <- polls_us_election_2016 %>% 
  filter(pollster == "Ipsos" & state == "U.S.") %>% 
  mutate(spread = rawpoll_clinton/100 - rawpoll_trump/100)
```

Without pollster effects (from using multiple pollsters), perhaps the theoretical and data-derived standard error is the same.

``` r
se <- one_pollster %>% 
  summarize(empirical = sd(spread), 
            theoretical = 2*sqrt(mean(spread)*(1-mean(spread))/min(samplesize)))
se
```

    ##    empirical theoretical
    ## 1 0.04025194  0.03256719

In fact, the actual standard error is more than we'd expect. A look at the histogram of the spreads reveals that the distribution isn't normal.

``` r
one_pollster %>% ggplot(aes(spread)) + 
  geom_histogram(binwidth = 0.01, color = "black")
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-49-1.png)

The following plots show us that the assumption we've made that *p* is fixed does not account for fluctuations that occur over time due to such factors as party conventions and revelations about private deeds.

``` r
polls_us_election_2016 %>%
  filter(state == "U.S." & enddate>="2016-07-01") %>%
  group_by(pollster) %>%
  filter(n()>=10) %>%
  ungroup() %>%
  mutate(spread = rawpoll_clinton/100 - rawpoll_trump/100) %>%
  ggplot(aes(enddate, spread)) + 
  geom_smooth(method = "loess", span = 0.1) + 
  geom_point(aes(color=pollster), show.legend = FALSE, alpha=0.6) 
```

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-50-1.png)

This variation implies that we need to include a time effect term to our model:

*Y*<sub>*i**j**t*</sub> = *d* + *b* + *h*<sub>*i*</sub> + +*b*<sub>*t*</sub> + *ϵ*<sub>*i*, *j*</sub>
 Pollsters also try to estimate trends, *f*(*t*), and include them in their predictions:

*Y*<sub>*i**j**t*</sub> = *d* + *b* + *h*<sub>*i*</sub> + +*b*<sub>*t*</sub> + *f*(*t*)+*ϵ*<sub>*i*, *j*</sub>
 The trends, *f*(*t*), are usually estimated for the proportion rather than the spread.

``` r
polls_us_election_2016 %>%
  filter(state == "U.S." & enddate>="2016-07-01") %>%
  select(enddate, pollster, rawpoll_clinton, rawpoll_trump) %>%
  rename(Clinton = rawpoll_clinton, Trump = rawpoll_trump) %>%
  gather(candidate, percentage, -enddate, -pollster) %>% 
  mutate(candidate = factor(candidate, levels = c("Trump","Clinton")))%>%
  group_by(pollster) %>%
  filter(n()>=10) %>%
  ungroup() %>%
  ggplot(aes(enddate, percentage, color = candidate)) +  
  geom_point(show.legend = FALSE, alpha=0.4)  + 
  geom_smooth(method = "loess", span = 0.15) +
  scale_y_continuous(limits = c(30,50))
```

    ## Warning: Removed 22 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 22 rows containing missing values (geom_point).

![](R_Inference_and_Modeling_Course_files/figure-markdown_github/unnamed-chunk-51-1.png)

### The t-Distribution

If the population data is known to be normal, we can use the t-distribution to estimate the increase in our standard error when our standard deviation is unknown (e.g. N is low). When using the t-distribution, we are estimating sigma, which introduces additional variability into our result.

We've used the following formula to represent the normal variable, *Z*.

$$Z = \\frac{\\bar{X} - d}{\\sigma/\\sqrt{N}}$$
 In practice, we have to estimate *σ*, which show by replacing it with *s*, which is estimated from the data.

$$Z = \\frac{\\bar{X} - d}{s/\\sqrt{N}}$$
 The theory tells us that as long as the larger population is normally distributed, we can use the t-distribution to describe *Z* with our estimated error *s*.

In the following example, we first show that the 95% CI calculated using a normal distribution captures *μ* less than 95% of the time when *N* = 15.

``` r
data(heights)
x <- heights %>% filter(sex == "Male") %>%
  .$height
set.seed(1)
mu <- mean(x)
N <- 15
B <- 10
res <- replicate(B, {
  i <- sample(x, N, replace = TRUE)
  interval <- c(qnorm(0.025, mean = mean(i), sd = sd(i)/sqrt(N)), qnorm(0.975, mean = mean(i), sd = sd(i)/sqrt(N)))
  between(mu, interval[1], interval[2])
})
mean(res)
```

    ## [1] 0.9

If we use the t-distribution instead, our 95% CI is more accurate, actually encompassing *μ* about 95% of the time.

``` r
mu <- mean(x)
set.seed(1)
N <- 15
B <- 10000
res <- replicate(B, {
  i <- sample(x, N, replace = TRUE)
  interval <- c(mean(i) - qt(0.975,df=N-1)*sd(i)/sqrt(N), mean(i) + qt(0.975,df=N-1)*sd(i)/sqrt(N))
  between(mu, interval[1], interval[2])
})
mean(res)
```

    ## [1] 0.9523

As you can see, the t-distribution accounts for the increase in variability that occurs with low N. Note the difference in how `qt()` functions, using only a proportion and degrees of freedom as its inputs. Multiplying by the standard error and adding or substracting from the sample mean produces the 95% and 5% boundaries of the CI, respectively.

> The t-distribution approaches normality as N increases. A t-distribution with N ≥ 30 looks essentially normal.

Look out for situations where the t-distribution is appropriate, such as when you are estimating *σ* and/or *N* is low.

Association Tests
-----------------