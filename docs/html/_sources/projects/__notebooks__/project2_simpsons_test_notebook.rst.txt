.. _proj2_simpsons_test:

A simple notebook for validating Simpson’s rule for a simple polynomial
-----------------------------------------------------------------------

First we run \_test_numerical_integration() to integrate the following
polynomial using both trapezoid and simpsons rule’s for integration.

:math:`f(x) = x^3 +5x^2 -2.125x - 1.1521`

And the exact result is:

:math:`\int_{-10}^{10}f(x)dx = 3310.291333333333`

Since Simpson’s rule integrates polynomials of up to and including
degree three exactly, then we should get exact results when using
Simpson’s rule here

Simpsons rule requires an even number of ‘bins’ but can be used with an
odd number of bins which we have also implemented - when using an odd
number of bins integration will be ‘close’ but not exact - whereas an
even number of bins will provide the exact result

See: Composite Simpson’s rule for irregularly spaced data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

https://en.wikipedia.org/wiki/Simpson%27s_rule

.. code:: ipython3

    # First we test the implemented numerical integration
    # functions to prove to ourselves that our methods
    # are not gibberish since latest code update
    from transportcorrection import _test_numerical_integration
    _test_numerical_integration()


.. parsed-literal::

    Now running numerical integration for various settings ...
    Simpsons rule with even N should be exact since it is of low enough polynomial order.
    Exact: 3310.291333333333
    Simpsons Rule Odd N: 3316.387964944368 | time: 0.00018858909606933594
    Simpsons Rule Even N: 3310.291333333333 | time: 6.198883056640625e-05
    Trapezoid Rule N=10: 3376.9579999999996 | time: 5.3882598876953125e-05
    Trapezoid Rule N=20: 3326.958 | time: 5.078315734863281e-05
    Trapezoid Rule N=30: 3317.6987407407405 | time: 5.14984130859375e-05
    Trapezoid Rule N=40: 3314.4579999999996 | time: 5.078315734863281e-05
    Trapezoid Rule N=50: 3312.9579999999996 | time: 5.1975250244140625e-05
    Trapezoid Rule N=1000: 3310.2980133533592 | time: 0.0002682209014892578
    Trapezoid Rule N=3000: 3310.292074568148 | time: 0.00017452239990234375
    Simpsons Rule N=3000: 3310.2913333333336 | time: 0.0015559196472167969
    Simpsons Rule N=2999: 3310.2913333338333 | time: 0.0015690326690673828

