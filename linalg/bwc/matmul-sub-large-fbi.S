#if 0
static inline void fill_buckets_indirect(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abcopy(x, sb[*q], z, 1);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}
#endif

        .text
	.p2align 4,,15
.globl fill_buckets_indirect
	.type	fill_buckets_indirect, @function
fill_buckets_indirect:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12

        movq %rsi,%rdi  // sb
        movq %rdx,%rsi  // z
        // rcx == 1
        // r8 == n
	testl	%r8d, %r8d
	je	.Ldone
.Lpre:
        movl %r8d, %eax
        andl $3, %eax
        je .Laligned
	movzbl	(%rcx), %eax
	leaq	(%rsi,%rax,8), %rsi
	movzbl	1(%rcx), %eax
	movq	(%rdi,%rax,8), %rdx
	movq	(%rsi), %rax
	movq	%rax, (%rdx)
	movzbl	1(%rcx), %eax
	addq	$2, %rcx
	addq	$8, (%rdi,%rax,8)
        subl    $1, %r8d
        ja .Lpre
.Laligned:
	testl	%r8d, %r8d
        je .Ldone
	.p2align 4,,10
	.p2align 3
        // store end pointer in rbx
        leaq    (%rcx,%r8,2), %rbx

/////////////////////////////////////////////////////////

.Louter_begin:
        // number of times the inner loop is run.
        // the inner loop scans 4 coefficients
        movq    $0x4, %rax
        leaq    (%rcx,%rax,8), %rbp
        movq    %rbx, %rdx
        cmpq    %rbp, %rbx
        // rbx rdx far fence
        // rbp close fence
        ja .Lshort
        movq    %rbp, %rdx
        // rbx far fence
        // rbp rdx close fence
.Lshort:
        movq    %rdx, %rbp

        xorq    %r9, %r9
.Lprefetch:
        leaq    (%rsi,%r9,0x1), %rdx
        prefetchnta      0x400(%rdx)
        addq    $0x40, %r9
        cmpq    $0x2000, %r9
        ja      .Lprefetch

        // the critical stride for our 256set, 16-way associative L2 data cache
        // is 0x10 * 0x40 = 0x400
        // for an average dj value of 2, we expect to read, in order:
        // 0x10(%rsi)
        // 0x20(%rsi)
        // ...
        // 0x100(%rsi)

.Linner_begin:
	movzbl	0(%rcx), %r8d
	movzbl	1(%rcx), %r9d
	movzbl	2(%rcx), %r10d
	movzbl	3(%rcx), %r11d
	movzbl	4(%rcx), %r12d
	movzbl	5(%rcx), %r13d
	movzbl	6(%rcx), %r14d
	movzbl	7(%rcx), %r15d

        addq    %r8, %r10
        addq    %r10, %r12
        addq    %r12, %r14

	leaq	(%rsi,%r14,8), %r14

	movq	(%rsi,%r8,8), %r8
	movq	(%rsi,%r10,8), %r10
	movq	(%rsi,%r12,8), %r12
        movq    %r14, %rsi
	movq	(%r14), %r14

	movq	(%rdi,%r9,8), %rdx
	movq	%r8, (%rdx)
	addq	$8, (%rdi,%r9,8)

	movq	(%rdi,%r11,8), %rax
	movq	%r10, (%rax)
	addq	$8, (%rdi,%r11,8)

	movq	(%rdi,%r13,8), %rdx
	movq	%r12, (%rdx)
	addq	$8, (%rdi,%r13,8)

	movq	(%rdi,%r15,8), %rax
	movq	%r14, (%rax)
	addq	$8, (%rdi,%r15,8)

	addq	$8, %rcx

        cmpq    %rcx,%rbp
	ja	.Linner_begin
.Linner_end:

        cmpq    %rcx,%rbx
	ja	.Louter_begin
.Louter_end:

/////////////////////////////////////////////////////////
.Ldone:
        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
    
	.size	fill_buckets_indirect, .-fill_buckets_indirect
	.p2align 4,,15
