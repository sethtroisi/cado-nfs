	.text
.Ltext0:
	.p2align 4,,15
.globl matmul_sliced_asm
	.type	matmul_sliced_asm, @function
matmul_sliced_asm:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12
        pushq %r11
        pushq %r10
        pushq %r9
        pushq %r8
        movq    %rsi, %rdi
	movzwl	2(%rcx), %eax
	movq	%rdx, %rsi
	movzwl	(%rcx), %edx
	addq	$4, %rcx
	salq	$16, %rax
	movq	%rax, %r9
	orq	%rdx, %r9
        leaq    (%rcx,%r9,4),%rbx
	je	.L2

        movq    $4, %rax
        subq    %r9, %rax
        andq    $3, %rax
        je .L30
        subq    $4, %rcx
        decq    %rax
        je .L31
        subq    $4, %rcx
        decq    %rax
        je .L32
        subq    $4, %rcx
        jmp .L33
.L3:
	movzwq	(%rcx), %r10
	movzwq	2(%rcx), %r11
	leaq	(%rsi,%r10,8), %rsi
	movq	(%rsi), %rdx
	xorq	%rdx, (%rdi,%r11,8)
.L31:
	movzwq	4(%rcx), %r10
	movzwq	6(%rcx), %r11
	leaq	(%rsi,%r10,8), %rsi
	movq	(%rsi), %rdx
	xorq	%rdx, (%rdi,%r11,8)
.L32:
	movzwq	8(%rcx), %r10
	movzwq	10(%rcx), %r11
	leaq	(%rsi,%r10,8), %rsi
	movq	(%rsi), %rdx
	xorq	%rdx, (%rdi,%r11,8)
.L33:
	movzwq	12(%rcx), %r10
	movzwq	14(%rcx), %r11
	leaq	(%rsi,%r10,8), %rsi
	movq	(%rsi), %rdx
	xorq	%rdx, (%rdi,%r11,8)

	addq	$16, %rcx
	cmpq	%rcx, %rbx
	je	.L2

.L30:
.L4:
	movzwq	(%rcx), %r8
	movzwq	2(%rcx), %r9
	movzwq	4(%rcx), %r10
	movzwq	6(%rcx), %r11
	movzwq	8(%rcx), %r12
	movzwq	10(%rcx), %r13
	movzwq	12(%rcx), %r14
	movzwq	14(%rcx), %r15
        salq    $3, %r8
        salq    $3, %r10
        salq    $3, %r12
        salq    $3, %r14


	movq	(%rsi,%r8), %rdx
	xorq	%rdx,(%rdi,%r9,8)
        addq    %r10, %r8
	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r11,8)
        addq    %r12, %r8
	movq	(%rsi,%r8), %rdx
	xorq	%rdx, (%rdi,%r13,8)
        addq    %r14, %r8
	movq	(%rsi,%r8), %rax
	xorq	%rax, (%rdi,%r15,8)
        addq    %r8, %rsi

	addq	$16, %rcx
	cmpq	%rcx, %rbx
        
	ja	.L4

.L2:
	movq	%rbx, %rax
        popq %r8
        popq %r9
        popq %r10
        popq %r11
        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
.LFE52:
	.size	matmul_sliced_asm, .-matmul_sliced_asm
