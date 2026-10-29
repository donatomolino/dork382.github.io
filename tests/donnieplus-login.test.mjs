import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Window } from 'happy-dom';
import { readFileSync, readdirSync } from 'node:fs';
import { dpReadCredentials } from '../public/donnieplus/donnieplus-auth.js';

test('password-manager values submit even without React input/change events', () => {
  const { document } = new Window();
  const form = document.createElement('form');
  form.innerHTML = '<input id="email" name="email"><input id="password" name="password" type="password">';
  form.querySelector('#email').value = '  qa@example.test  ';
  form.querySelector('#password').value = ' fixture-only password ';
  // Controlled state deliberately remains empty: this is the autofill case.
  assert.deepEqual(dpReadCredentials(form, '', ''), { email: 'qa@example.test', password: ' fixture-only password ' });
  form.querySelector('#password').value = '';
  assert.equal(dpReadCredentials(form, 'old@example.test', 'old password').password, '', 'cleared DOM fields must not reuse old credentials');
});

test('modern and legacy login bundles provide native submission and DOM-value synchronization', () => {
  const base = new URL('../public/donnieplus/assets/', import.meta.url);
  for (const name of readdirSync(base).filter(name => /^Account-/.test(name))) {
    const source = readFileSync(new URL(name, base), 'utf8');
    assert.match(source, /dpReadCredentials\([IL]\.current,b,[CS]\)/);
    assert.ok(source.includes('id:e,name:e,className:`field__input'));
    assert.match(source, /`button`,\{type:`submit`,className:`btn btn--primary auth__submit/);
    assert.ok(source.includes('disabled:'));
    assert.ok(source.includes('onFocusCapture:()=>{const dpFields=dpReadCredentials'));
    assert.ok(source.includes('onBlurCapture:()=>{const dpFields=dpReadCredentials'));
  }
});
