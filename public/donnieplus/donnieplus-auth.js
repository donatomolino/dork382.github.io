// Read at submission time: password managers may update the DOM without
// delivering the same input/change events as manual typing.
export function dpReadCredentials(form, email, password) {
  const emailField = form && form.querySelector('input[id="email"]');
  const passwordField = form && form.querySelector('input[id="password"]');
  return { email: String(emailField ? emailField.value : email || '').trim(), password: String(passwordField ? passwordField.value : password || '') };
}
